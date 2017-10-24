###
### Simulation for survival probability
### Alison Ketz 10/11/2017
###

#preliminaries

rm(list=ls())

#functions

logit = function(x){log(x/(1-x))}
expit = function(x){exp(x)/(1+exp(x))}

#setseed
set.seed(121915)

#libraries

library(rjags)
library(R2jags)
library(coda)
library(dclone)
library(parallel)
library(foreach)
library(doParallel)
library(ggplot2)


############################################################################################################################################################
###
### Multistate Mark Recapture with living disease state + dead disease state
###
############################################################################################################################################################
# preliminaries
# 
# rm(list=ls())
# 
# #functions
# 
# logit = function(x){log(x/(1-x))}
# expit = function(x){exp(x)/(1+exp(x))}
# 
# #setseed
# set.seed(121915)
# 
# #libraries
# 
# library(rjags)
# library(R2jags)
# library(coda)
# 


# Specify model in JAGS language
sink("ms.disease.jags.R")
cat("
    model {
      
      # -------------------------------------------------
      # Parameters
      # s1: survival probability neg
      # s2: survival probability pos
      # psi: transition probability from neg to pos
      # p1: recapture probability for neg
      # p2: recapture probability for pos
      # -------------------------------------------------
      
      # -------------------------------------------------
      # States (S):
      # 1 alive neg
      # 2 alive pos
      # 3 dead neg
      # 4 dead pos
      # Observations (O):  
      # 1 seen neg 
      # 2 seen pos
      # 3 dead recovered neg
      # 4 dead recovered pos
      # 3 not seen
      # -------------------------------------------------
        
      # Priors and constraints
      for (t in 1:(n.occasions-1)){
        p1[t] ~ dunif(0,1)
        p2[t] ~ dunif(0,1)
        p3[t] ~ dunif(0,1)
        p4[t] ~ dunif(0,1)
      }
      
      s1 ~ dbeta(10,1) #annual survival cwd-
      s2 ~ dbeta(1,10) #annual survival cwd+
      psi ~ dunif(0, 1)    # Priors for disease state transition
      
      # Define state-transition and observation matrices
      for (i in 1:nind){  
        # Define probabilities of state S(t+1) given S(t)
        for (t in f[i]:(n.occasions-1)){
          ps[1,i,t,1] <- s1 * (1-psi)
          ps[1,i,t,2] <- s1 * psi
          ps[1,i,t,3] <- (1-s1)*(1-psi)
          ps[1,i,t,4] <- (1-s1)*psi
          
          ps[2,i,t,1] <- 0 
          ps[2,i,t,2] <- s2
          ps[2,i,t,3] <- 0
          ps[2,i,t,4] <- 1-s2
          
          ps[3,i,t,1] <- 0
          ps[3,i,t,2] <- 0
          ps[3,i,t,3] <- 1
          ps[3,i,t,4] <- 0
          
          ps[4,i,t,1] <- 0
          ps[4,i,t,2] <- 0
          ps[4,i,t,3] <- 0
          ps[4,i,t,4] <- 1
          
          # Define probabilities of O(t) given S(t)
          po[1,i,t,1] <- p1[t]
          po[1,i,t,2] <- 0
          po[1,i,t,3] <- 0
          po[1,i,t,4] <- 0
          po[1,i,t,5] <- 1-p1[t]
          
          po[2,i,t,1] <- 0
          po[2,i,t,2] <- p2[t]
          po[2,i,t,3] <- 0
          po[2,i,t,4] <- 0
          po[2,i,t,5] <- 1-p2[t]
          
          po[3,i,t,1] <- 0
          po[3,i,t,2] <- 0
          po[3,i,t,3] <- p3[t]
          po[3,i,t,4] <- 0
          po[3,i,t,5] <- 1-p3[t]
          
          po[4,i,t,1] <- 0
          po[4,i,t,2] <- 0
          po[4,i,t,3] <- 0
          po[4,i,t,4] <- p4[t]
          po[4,i,t,5] <- 1-p4[t]
        } #t
      } #i
  
      # Likelihood 

      for (i in 1:nind){
        # Define latent state at first capture
        z[i,f[i]] <- y[i,f[i]]
        for (t in (f[i]+1):n.occasions){
          # State process: draw S(t) given S(t-1)
          z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
          # Observation process: draw O(t) given S(t)
          y[i,t] ~ dcat(po[z[i,t], i, t-1,])
        } #t
      } #i
    }
    ",fill = TRUE)
sink()

# Write function to simulate data


# Define function to simulate multistate capture-recapture data
simul.ms <- function(PSI.STATE, PSI.OBS, marked, unobservable = NA){
  
  # Unobservable: number of state that is unobservable
  n.occasions <- dim(PSI.STATE)[4] + 1
  
  CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
  
  # Define a vector with the occasion of marking
  mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked))
  g <- colSums(marked)
  for (s in 1:dim(PSI.STATE)[1]){
    if (g[s]==0) next  # To avoid error message if nothing to replace
    mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <-
      rep(1:n.occasions, marked[1:n.occasions,s])
  } #s
  for (i in 1:sum(marked)){
    for (s in 1:dim(PSI.STATE)[1]){
      if (mark.occ[i,s]==0) next
      first <- mark.occ[i,s]
      CH[i,first] <- s
      CH.TRUE[i,first] <- s
    } #s
    for (t in (first+1):n.occasions){
      # Multinomial trials for state transitions
      if (first==n.occasions) next
      state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
      CH.TRUE[i,t] <- state
      # Multinomial trials for observation process
      event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
      CH[i,t] <- event
    } #t
  } #i
  # Replace the NA and the highest state number (dead) in the file by 0
  CH[is.na(CH)] <- 0
  CH[CH==dim(PSI.STATE)[1]] <- 0
  CH[CH==unobservable] <- 0
  id <- numeric(0)
  for (i in 1:dim(CH)[1]){
    z <- min(which(CH[i,]!=0))
    ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id))
  }
  return(list(CH=CH[-id,], CH.TRUE=CH.TRUE[-id,]))
  # CH: capture histories to be used
  # CH.TRUE: capture histories with perfect observation
}


# Function to create known latent states z
known.state.ms <- function(ms, notseen){
  # notseen: label for ?not seen?
  state <- ms
  state[state==notseen] <- NA
  for (i in 1:dim(ms)[1]){
    m <- min(which(!is.na(state[i,])))
    state[i,m] <- NA
  }
  return(state)
}


# Function to create initial values for unknown z
ms.init.z <- function(ch, f, z.known,n.occasions){
  for(i in 1:dim(ch)[1]){
    for(t in 1:n.occasions){
      if(ch[i,t]==5)ch[i,t]=NA
    }
  }
  
  for(i in 1:dim(ch)[1]){
    for(t in f[i]:n.occasions){
      if(is.na(ch[i,t]))ch[i,t]=ch[i,t-1]
    }
  } 
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  for(i in 1:dim(ch)[1]){
    for(t in 1:n.occasions){
      if(!is.na(z.known[i,t]))ch[i,t]=NA
    }
  }  
  return(ch)
}


###
### Define overall simulation function can change values of p1.true
### returns draws from s1
###
ni=50000 #set number of mcmc iterations
nc = 3

simulation.survive=function(p1.true=.01,ni=50000,nc=3){
  n.occasions = 5
  n.states = 4
  n.obs = 5
  s1= .6
  s2 = .3
  psi = .1 #transition probability from susceptible (neg) to infected (pos)
  p1 = rep(p1.true,n.occasions-1)
  p2 = rep(.01,n.occasions-1)
  p3 = c(.25,0.25,.9,.9)
  p4 = c(.25,.25,.9,.9)
  marked = matrix(NA,ncol=n.states,nrow = n.occasions)
  marked[,1] = rep(100,n.occasions) #marked negative,
  marked[,2] = rep(30,n.occasions) #marked positive
  marked[,3] = rep(0,n.occasions) #marked dead negative
  marked[,4] = rep(0,n.occasions) #marked dead negative
  
  
  # define matrices with survival, transition, and recapture probabilites
  # this is a 4 dimensional matrix
  # dimension 1: state of 'departure'
  # dimension 2: state of 'arrival'
  # dimension 3: individual
  # dimension 4: time
  
  #state process matrix
  totrel = sum(marked)*(n.occasions-1)
  PSI.STATE = array(NA,dim=c(n.states,n.states,totrel,n.occasions-1))
  for(i in 1:totrel){
    for(t in 1:(n.occasions-1)){
      PSI.STATE[,,i,t] = matrix(c(s1*(1-psi),s1*psi,(1-s1)*(1-psi),(1-s1)*psi,
                                  0,s2,0,1-s2,
                                  0,0,1,0,
                                  0,0,0,1),nrow=n.states,byrow=TRUE)
    }
  }
  
  #observation process matrix
  PSI.OBS = array(NA,dim=c(n.states,n.obs,totrel,n.occasions-1))
  for(i in 1:totrel){
    for(t in 1:(n.occasions-1)){
      PSI.OBS[,,i,t] = matrix(c(p1[t],0,0,0,1-p1[t],
                                0,p2[t],0,0,1-p2[t],
                                0,0,p3[t],0,1-p3[t],
                                0,0,0,p4[t],1-p4[t]),nrow=n.states,byrow=TRUE)
    }
  }
  
  
  #execute simulation
  
  ms.dis.sim =  simul.ms(PSI.STATE,PSI.OBS,marked)
  CH = ms.dis.sim$CH
  
  #calculate prevalence to make sure probability of transmission isn't off-base
  # 
  # alive.pos=which(ms.dis.sim$CH.TRUE==2,TRUE)
  # alive.neg=which(ms.dis.sim$CH.TRUE==1,TRUE)
  # prev.vec=c()
  # for(t in 1:n.occasions){
  #   prev.vec=c(prev.vec,sum(alive.pos[,2]==t)/(sum(alive.pos[,2]==t)+sum(alive.neg[,2]==t)))
  # }
  # prev.vec
  
  #compute vector with occasion of first capture
  get.first=function(x)min(which(x!=0))
  f = apply(CH,1,get.first)
  
  # 1 = seen alive and neg, 2 = seen alive and pos, 3,dead and neg, 4 dead and pos, 5 = not seen
  rCH <- CH          # Recoded CH
  rCH[rCH==0] <- 5
  
  #generating initial values
  z.known = known.state.ms(rCH, 5)
  z.init = ms.init.z(rCH, f,z.known,n.occasions)

  # Bundle data
  jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = z.known)
  
  
  # Initial values
# 
#   inits1=list(s1=runif(1,.3,1),s2=runif(1,0,.4), psi = runif(1, 0, 1), p1=p1+.01,p2=p2+.01,p3=p3+.01,p4=p4+.01, z = z.init)
#   inits2=list(s1=runif(1,.3,1),s2=runif(1,0,.4), psi = runif(1, 0, 1), p1=p1,p2=p2,p3=p3,p4=p4, z = z.init)
#   inits3=list(s1=runif(1,.3,1),s2=runif(1,0,.4), psi = runif(1, 0, 1), p1=p1-.01,p2=p2-.01,p3=p3-.01,p4=p4-.01, z = z.init)
#   
  inits <- list(list(s1=runif(1,.3,1),s2=runif(1,0,.4), psi = runif(1, 0, 1), p1=rep(.1,4),p2=rep(.1,4),p3=rep(.1,4),p4=rep(.1,4), z = z.init),
                             list(s1=runif(1,.3,1),s2=runif(1,0,.4), psi = runif(1, 0, 1), p1=rep(.3,4),p2=rep(.3,4),p3=rep(.3,4),p4=rep(.3,4), z = z.init),
                             list(s1=runif(1,.3,1),s2=runif(1,0,.4), psi = runif(1, 0, 1), p1=rep(.2,4),p2=rep(.2,4),p3=rep(.2,4),p4=rep(.2,4), z = z.init))

  
  # Parameters monitore4
  parameters <- c("s1","s2" ,"psi", "p1","p2","p3","p4")
  
  # MCMC settings
  nt = 1
  nb = 10000

  # call parallel version of jags using dclone
  cl=makeCluster(3)
  ms.cwd=jags.parfit(cl, data=jags.data, params=parameters, model="ms.disease.jags.R",inits=inits,n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
  stopCluster(cl)
  
  # Call JAGS (Run Time 6 min) using R2jags
  # ms.cwd <- jags(jags.data, inits, parameters, "ms.disease.jags.R", n.chains = 3, n.iter = 10000)
  gel= gelman.diag(ms.cwd,multivariate = FALSE)

  ms.cwd.iter=as.data.frame(rbind(ms.cwd[[1]],ms.cwd[[2]],ms.cwd[[3]]))

  return(list(s1.iter=ms.cwd.iter$s1,psi.iter=ms.cwd.iter$psi,gel=gel))
  
}#end simulation call

out=simulation.survive(.5)


p1.possible=c(0,.01,seq(0,1,by=.1)[-1])

s1.results=matrix(NA,nrow=ni*nc, ncol=length(p1.possible))
psi.results=matrix(NA,nrow=ni*nc, ncol=length(p1.possible))
gelman.results = vector(mode = "list", length = length(p1.possible))
out=list(s1.results,psi.results,gelman.results)

base=2
numcores=detectCores()-1
clmain=makeCluster(numcores)
registerDoParallel(clmain)

stopCluster(clmain)

# for(i in 1:p1.possible){
foreach(i = 1:2,.export="p1.possible"){
  
}
  results=simulation.survive(p1.possible[i],ni,nc)
  out[[1]][,i] = results[[1]]
  out[[2]][,i] = results[[2]]
  out[[3]][[i]] = results[[3]]
}
head(out[[1]])
head(out[[2]])
head(out[[3]])

head(out$psi.results)


# Call JAGS from R (BRT 8 min)
ms.cwd <- jags.(jags.data, inits, parameters, "ms.disease.jags.R", n.chains = 3, n.thin = 1, n.iter = 10000, n.burnin = 5000, n.cluster=3,working.directory = getwd())



ms.cwd.mcmc=as.mcmc(ms.cwd)
ms.cwd.iter=as.data.frame(rbind(ms.cwd.mcmc[[1]],ms.cwd.mcmc[[2]],ms.cwd.mcmc[[3]]))
names(ms.cwd.iter)=c("deviance","psi","s1","s2","p11","p12","p13","p14","p21","p22","p23","p24","p31","p32","p33","p34","p41","p42","p43","p44")

quantile(ms.cwd.iter$s1,.025)
hist(ms.cwd.iter$s1,main=expression(s[1]),xlab=expression(s[1]))
abline(v=s1,lwd=3,lty=2,col=2)
abline(v=quantile(ms.cwd.iter$s1,.025),lwd=3,lty=2,col=2)
abline(v=quantile(ms.cwd.iter$s1,.025),lwd=3,lty=2,col=2)
print(ms.cwd.iter)

psi
s1
s2
p1
p2
p3
p4




traceplot(ms.cwd)

ms.cwd.fit=as.mcmc(ms.cwd)
gelman.diag(ms.cwd.fit)

dim(ms.cwd.fit[[1]])

print(ms.cwd)



