

############################################################################################################################################################
###
### CJS simulation
###
############################################################################################################################################################


# Define parameter values
n.occasions = 5 # number of capture occasions
marked = rep(100, n.occasions-1) #annual number of newly marked individuals
mark.occ = rep(1:length(marked),marked[1:length(marked)])#define a vector with marking occasion

mark.occ


#define 2 survival probabilities
phi0=rep(0.45,n.occasions-1) #assume survival is .45 of healthy individuals
phi1=rep(0.3,n.occasions-1) #assume survival is .3 of cwd individuals

beta0 = logit(phi0)
beta1 = logit(phi1) - beta0

#probability of recapture
p=rep(.25,n.occasions-1) #assume recapture probability is .25

#disease status, 1 = cwd+
#beta0+beta1*x
psi=.3 #prevalence  = prob(cwd+)
delta = .2 #transmission prob

#Define matrix for Disease status, 1 = +, 0 = -
CWD = matrix (0,ncol=n.occasions, nrow=sum(marked))

for(i in 1:sum(marked)){
  CWD[i,mark.occ[i]] =  rbinom(1,1,psi)
  for(t in (mark.occ[i]+1):n.occasions){
    if(CWD[i,t-1] == 1){
      CWD[i,t]=1
    }
    else{
      CWD[i,t] = rbinom(1,1,delta)
    }
  }
}

#Define matrices with survival and recapture probabilities
PHI = matrix (phi0,ncol=n.occasions, nrow=sum(marked))

for(i in 1:sum(marked)){
  for (t in 1:n.occasions){
    if(CWD[i,t] == 1)PHI[i,t]=phi1
  }
}

P = matrix (p,ncol=n.occasions, nrow=sum(marked))

#Define function to simulate a capture history
simul.cjs = function(PHI,P,marked){
  n.occasions=dim(PHI)[2]+1
  CH = matrix(0,ncol=n.occasions,nrow=sum(marked))
  
  #define a vector with marking occasion
  mark.occ = rep(1:length(marked),marked[1:length(marked)])
  
  #Fill the CH matrix
  for (i in 1:sum(marked)){
    CH[i,mark.occ[i]] = 1
    if (mark.occ[i]==n.occasions) next
    for(t in (mark.occ[i]+1):n.occasions){
      
      #Bernoulli random trial - does individual survive occassion
      sur = rbinom(1,1,PHI[i,t-1])
      if(sur==0)break
      rp = rbinom(1,1,P[i,t-1])
      if(rp==1)CH[i,t] = 1
    }#t
  }#i
  return(CH)
}

CH = simul.cjs(PHI,P,marked)


#create vector with occassion of first marking
get.first = function(x)min(which(x!=0))
f = apply(CH,1,get.first)
#note: f = mark.occ


#specify model with jags code

sink("cjs.jags.R")
cat("
    model{
    #priors and constraints
    for(i in 1:nind){
    for(t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- beta0 + CWD[i,t]*beta1
    p[i,t] <- mean.p
    }#t
    }#i
    
    mean.p ~ dunif(0,1)
    beta0 ~ dnorm(0,1)
    beta1 ~ dnorm(0,1)
    
    #likelihood
    for(i in 1:nind){
    z[i,f[i]] <- 1 #define latent state presence at first capture
    for( t in (f[i]+1):n.occasions){
    #state process
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t] <- phi[i,t-1]*z[i,t-1]
    #observation process
    y[i,t] ~ dbern(mu2[i,t])
    mu2[i,t] <-p[i,t-1]*z[i,t]
    
    }#t
    
    }#i
    
    }#end model
    
    
    ",fill=TRUE)
sink()


#set data
data.jags = list(y=CH,f=f,nind=dim(CH)[1],n.occasions = dim(CH)[2],CWD = CWD)

#create a matrix of initial values for latent state z

known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}

inits=function(){list( z = known.state.cjs(CH),mean.p = p[1], beta0 = beta0[1]+.05,beta1 = beta1[1]+.05)}
parameters = c("mean.p","beta0","beta1")

#MCMC settings
ni=10000
nt=1
nb = 5000
nc = 3

model.fit = jags(data.jags, inits, parameters, "cjs.jags.R", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(model.fit,digits=3)

############################################################################################################################################################
###
### Mark-Recovery survival model - as state space model
###
############################################################################################################################################################

rm(list=ls())

logit = function(x){log(x/(1-x))}
expit = function(x){exp(x)/(1+exp(x))}

n.occasions = 4 #number of release occasions
marked = rep(100,n.occasions) #annual number of marked individuals
mark.occ = rep(1:n.occasions,marked)


s = rep(.5,n.occasions) # true survival probability
r = rep(.25,n.occasions) # recovery probability

S = matrix(s,ncol=n.occasions,nrow=sum(marked))
R = matrix(r,ncol=n.occasions,nrow=sum(marked))


#define function to simulate mark-recovery data

simul.mr = function(S,R,marked){
  n.occasions = dim(S)[2]
  MR = matrix(NA, ncol=n.occasions+1,nrow=sum(marked))
  
  #define vector of occassions marked
  mark.occ = rep(1:n.occasions,marked)
  
  #Fill the MR matrix
  for(i in 1:sum(marked)){
    MR[i,mark.occ[i]]=1 #release/collared occasion
    for(t in mark.occ[i]:n.occasions){
      sur = rbinom(1,1,S[i,t])
      if(sur==1)next #if still alive move to the next occasion
      rp = rbinom(1,1,R[i,t])
      if(rp==0){
        MR[i,t+1] = 0
        break
      }
      if(rp==1){
        MR[i,t+1]=1
        break
      }
    }#t
  }#i
  MR[which(is.na(MR))]=0
  return(MR)    
  
}


MR = simul.mr(S,R,marked)

#Analysis of model

get.first=function(x)min(which(x!=0))
f = apply(MR,1,get.first)

#Specify JAGS model
sink("mr.jags.R")
cat("
    model{
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    s[i,t] <- mean.s
    r[i,t] <- mean.r
    } #t
    } #i
    mean.s ~ dunif(0, 1)          # Prior for mean survival
    mean.r ~ dunif(0, 1)          # Prior for mean recapture
    
    #Likelihood
    for(i in 1:nind){
    #Define latent state at first capture
    z[i,f[i]] = 1
    for(t in (f[i]+1):n.occasions){
    #State process
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t] <- s[i,t-1]*z[i,t-1]
    
    #Observation process
    
    y[i,t] ~ dbern(mu2[i,t])
    mu2[i,t] <- r[i,t-1]*(z[i,t-1]-z[i,t])
    }#t
    }#i
    
    }#model
    
    
    ",fill=TRUE)
sink()


# Define function to create a matrix with information about known latent state z
known.state.mr <- function(mr){
  state <- matrix(NA, nrow = dim(mr)[1], ncol = dim(mr)[2])
  rec <- which(rowSums(mr)==2)
  for (i in 1:length(rec)){
    n1 <- min(which(mr[rec[i],]==1))
    n2 <- max(which(mr[rec[i],]==1))
    state[rec[i],n1:n2] <- 1
    state[rec[i],n1] <- NA
    state[rec[i],n2:dim(mr)[2]] <- 0
  }
  return(state)
}

# Bundle data
data.jags <- list(y = MR, f = f, nind = dim(MR)[1], n.occasions = dim(MR)[2], z = known.state.mr(MR))

# Define function to create a matrix of initial values for latent state z
mr.init.z <- function(mr){
  ch <- matrix(NA, nrow = dim(mr)[1], ncol = dim(mr)[2])
  rec <- which(rowSums(mr)==1)
  for (i in 1:length(rec)){
    n1 <- which(mr[rec[i],]==1)
    ch[rec[i],n1:dim(mr)[2]] <- 0
    ch[rec[i],n1] <- NA
  }
  return(ch)
}

# Initial values
inits <- function(){list(z = mr.init.z(MR), mean.s = runif(1, 0, 1), mean.r = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.s", "mean.r")

# MCMC settings
ni <- 10000
nt <- 1
nb <- 5000
nc <- 3

mr.ss.fit <- jags(data.jags, inits, parameters, "mr.jags.R", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(mr.ss.fit, digits = 3)

s
r

############################################################################################################################################################
###
### Mark-Recovery survival model - As multinomial
###
############################################################################################################################################################


# 8.3. The mark-recovery model fitted with the multinomial likelihood 

# Define function to create an m-array based for mark-recovery (MR) data
marray.dead <- function(MR){
  nind <- dim(MR)[1]
  n.occasions <- dim(MR)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Create vector with occasion of marking 
  get.first <- function(x) min(which(x!=0))
  f <- apply(MR, 1, get.first)
  # Calculate the number of released individuals at each time period
  first <- as.numeric(table(f))
  for (t in 1:n.occasions){
    m.array[t,1] <- first[t]
  }
  # Fill m-array with recovered individuals
  rec.ind <- which(apply(MR, 1, sum)==2)
  rec <- numeric()
  for (i in 1:length(rec.ind)){
    d <- which(MR[rec.ind[i],(f[rec.ind[i]]+1):n.occasions]==1)
    rec[i] <- d + f[rec.ind[i]]
    m.array[f[rec.ind[i]],rec[i]] <- m.array[f[rec.ind[i]],rec[i]] + 1
  }
  # Calculate the number of individuals that are never recovered
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1]-sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

marr <- marray.dead(MR)

# Specify model in BUGS language
sink("mr.mnl.jags.R")
cat("
    model {
    
    # Priors and constraints
    for (t in 1:n.occasions){
    s[t] <- mean.s
    r[t] <- mean.r
    }
    mean.s ~ dunif(0, 1)              # Prior for mean survival
    mean.r ~ dunif(0, 1)              # Prior for mean recovery
    
    # Define the multinomial likelihood
    for (t in 1:n.occasions){
    marr[t,1:(n.occasions+1)] ~ dmulti(pr[t,], rel[t])
    }
    
    # Define the cell probabilities of the m-array
    # Main diagonal
    for (t in 1:n.occasions){
    pr[t,t] <- (1-s[t])*r[t]
    
    # Above main diagonal
    for (j in (t+1):n.occasions){
    pr[t,j] <- prod(s[t:(j-1)])*(1-s[j])*r[j]
    } #j
    
    # Below main diagonal
    for (j in 1:(t-1)){
    pr[t,j] <- 0
    } #j
    } #t
    
    # Last column: probability of non-recovery
    for (t in 1:n.occasions){
    pr[t,n.occasions+1] <- 1-sum(pr[t,1:n.occasions])
    } #t
    } #model
    ",fill = TRUE)
sink()

# Bundle data
jags.data <- list(marr = marr, n.occasions = dim(marr)[2]-1, rel = rowSums(marr))

# Initial values
inits <- function(){list(mean.s = runif(1, 0, 1), mean.r = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.s", "mean.r")

# MCMC settings
ni <- 10000
nt <- 1
nb <- 5000
nc <- 3

# Call JAGS from R (BRT <1 min)
mr <- jags(jags.data, inits, parameters, "mr.mnl.jags.R", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(mr, digits = 3)




############################################################################################################################################################
###
### Mark-Recovery survival model - As multinomial - with age
###
############################################################################################################################################################

# 8.3.2. Age-dependent parameters
n.occasions <- 5                   # Number of occasions
marked.j <- rep(200, n.occasions)   # Annual number of newly marked young
marked.a <- rep(20, n.occasions)    # Annual number of newly marked adults
sjuv <- 0.3                         # Juvenile survival probability
sad <- 0.8                          # Adult survival probability
rjuv <- 0.25                        # Juvenile recovery probability
rad <- 0.15                         # Adult recovery probability
sj <- c(sjuv, rep(sad, n.occasions-1))
rj <- c(rjuv, rep(rad, n.occasions-1))

# Define matrices with survival and recovery probabilities
SJ <- matrix(0, ncol = n.occasions, nrow = sum(marked.j))
for (i in 1:length(marked.j)){
  SJ[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:n.occasions] <- matrix(rep(sj[1:(n.occasions-i+1)],marked.j[i]), ncol = n.occasions-i+1, byrow = TRUE)
}
SA <- matrix(sad, ncol = n.occasions, nrow = sum(marked.a))
RJ <- matrix(0, ncol = n.occasions, nrow = sum(marked.j))
for (i in 1:length(marked.j)){
  RJ[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:n.occasions] <- matrix(rep(rj[1:(n.occasions-i+1)],marked.j[i]), ncol = n.occasions-i+1, byrow = TRUE)
}
RA <- matrix(rad, ncol = n.occasions, nrow = sum(marked.a))

# Execute simulation function
MRj <- simul.mr(SJ, RJ, marked.j)
MRa <- simul.mr(SA, RA, marked.a)

# Summarize data in m-arrays
marr.j <- marray.dead(MRj)
marr.a <- marray.dead(MRa)

# Specify model in BUGS language
sink("mr-mnl-age.jags")
cat("
    model {
    
    # Priors and constraints
    for (t in 1:n.occasions){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    rj[t] <- mean.rj
    ra[t] <- mean.ra
    }
    mean.sj ~ dunif(0, 1)              # Prior for mean juv. survival
    mean.sa ~ dunif(0, 1)              # Prior for mean ad. survival
    mean.rj ~ dunif(0, 1)              # Prior for mean juv. recovery
    mean.ra ~ dunif(0, 1)              # Prior for mean ad. recovery
    
    # Define the multinomial likelihoods
    # Calculate the number of birds released each year
    for (t in 1:n.occasions){
    marr.j[t,1:(n.occasions+1)] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:(n.occasions+1)] ~ dmulti(pr.a[t,], rel.a[t])
    }
    # Define the cell probabilities of the juvenile m-array
    # Main diagonal
    for (t in 1:n.occasions){
    pr.j[t,t] <- (1-sj[t])*rj[t]
    # Further above main diagonal
    for (j in (t+2):n.occasions){
    pr.j[t,j] <- sj[t]*prod(sa[(t+1):(j-1)])*(1-sa[j])*ra[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.j[t,j] <- 0
    } #j
    } #t
    for (t in 1:(n.occasions-1)){
    # One above main diagonal
    pr.j[t,t+1] <- sj[t]*(1-sa[t+1])*ra[t+1] 
    } #t
    # Last column: probability of non-recovery
    for (t in 1:n.occasions){
    pr.j[t,n.occasions+1] <- 1-sum(pr.j[t,1:n.occasions])
    } #t
    # Define the cell probabilities of the adult m-array
    # Main diagonal
    for (t in 1:n.occasions){
    pr.a[t,t] <- (1-sa[t])*ra[t]
    # Above main diagonal
    for (j in (t+1):n.occasions){
    pr.a[t,j] <- prod(sa[t:(j-1)])*(1-sa[j])*ra[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recovery
    for (t in 1:n.occasions){
    pr.a[t,n.occasions+1] <- 1-sum(pr.a[t,1:n.occasions])
    } #t
    }
    ",fill = TRUE)
sink()

# Bundle data
jags.data <- list(marr.j = marr.j, marr.a = marr.a, n.occasions = dim(marr.j)[2]-1, rel.j = rowSums(marr.j), rel.a = rowSums(marr.a))

# Initial values
inits <- function(){list(mean.sj = runif(1, 0, 1), mean.sa = runif(1, 0, 1), mean.rj = runif(1, 0, 1), mean.ra = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.sj", "mean.rj", "mean.sa", "mean.ra")

# MCMC settings
ni <- 5000
nt <- 6
nb <- 2000
nc <- 3

# Call JAGS from R (BRT <1 min)
mr.age <- jags(jags.data, inits, parameters, "mr-mnl-age.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(mr.age, digits = 3)





##################

n.occasions = 5

marked.neg = rep(70,n.occasions)
marked.pos = rep(30,n.occasions)

s.neg = .6 #cwd- survival probability
s.pos = .3 #cwd+ survival probability

r.neg = .9 #cwd- recovery probability
r.pos = .9 #cwd+ recovery probability 

sneg=c(rep(s.neg,floor(n.occasions/2)),rep(spos,n.occasions-floor(n.occasions/2)))
spos = rep(s.pos,n.occasions)

# Define matrices with survival and recovery probabilities
SNEG <- matrix(0, ncol = n.occasions, nrow = sum(marked.neg))
for (i in 1:length(marked.neg)){
  SNEG[(sum(marked.neg[1:i])-marked.neg[i]+1):sum(marked.neg[1:i]),i:n.occasions] <- matrix(rep(sneg[1:(n.occasions-i+1)],marked.neg[i]), ncol = n.occasions-i+1, byrow = TRUE)
}
head(SNEG)
dim(SNEG)



# 8.3.2. Age-dependent parameters
n.occasions <- 5                   # Number of occasions
marked.j <- rep(200, n.occasions)   # Annual number of newly marked young
marked.a <- rep(20, n.occasions)    # Annual number of newly marked adults
sjuv <- 0.3                         # Juvenile survival probability
sad <- 0.8                          # Adult survival probability
rjuv <- 0.25                        # Juvenile recovery probability
rad <- 0.15                         # Adult recovery probability
sj <- c(sjuv, rep(sad, n.occasions-1))
rj <- c(rjuv, rep(rad, n.occasions-1))

# Define matrices with survival and recovery probabilities
SJ <- matrix(0, ncol = n.occasions, nrow = sum(marked.j))
for (i in 1:length(marked.j)){
  SJ[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:n.occasions] <- matrix(rep(sj[1:(n.occasions-i+1)],marked.j[i]), ncol = n.occasions-i+1, byrow = TRUE)
}
head(SJ)
dim(SJ)


SA <- matrix(sad, ncol = n.occasions, nrow = sum(marked.a))
RJ <- matrix(0, ncol = n.occasions, nrow = sum(marked.j))
for (i in 1:length(marked.j)){ 
  RJ[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:n.occasions] <- matrix(rep(rj[1:(n.occasions-i+1)],marked.j[i]), ncol = n.occasions-i+1, byrow = TRUE)
}
RA <- matrix(rad, ncol = n.occasions, nrow = sum(marked.a))
# Execute simulation function
MRj <- simul.mr(SJ, RJ, marked.j)
MRa <- simul.mr(SA, RA, marked.a)

# Summarize data in m-arrays
marr.j <- marray.dead(MRj)
marr.a <- marray.dead(MRa)

