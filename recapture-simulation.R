###
### IPM simulation for precision to determine if new collars better than retesting already collared individual
###

set.seed(Sys.time())

#number of surveys
n.occassions = 10

#number classes
n.class = 3

#Total population size across occassions
N.total = rbinom(n.occassions,20000,.9)
  
#harvest rates females (from Jenelle)
h.rate=.26

#proportion of each class
class.prop = c(.4,.4,.2)

#Process model
#population sizes by classes
N.vec = rmultinom(n.occassions,N.total,class.prop) #dim(N.vec) = [n.class,n.occasions]

#Transmission
#set inciduence
psi=.3


#Estimating survival using CJS model

#set survival
s1 = .25 #fawns-females
s2 = .83 #does - healthy
s3 = .41 #does - sick

n.states=n.class+1
n.obs = 3

#probability of recapture
p=.5

#number of marked individuals
marked = t(matrix(rep(c(80,100,20,0),n.occassions),nr=n.states))

#total number of marked
totalrel = sum(marked)*(n.occassions-1)

#survival state process matrix 
S.STATE = array(NA,dim=c(n.states,n.states,totalrel,n.occassions-1))

for (i in 1:totalrel){
  for(t in 1:(n.occassions-1)){
    S.STATE[,,i,t]=matrix(c(
      0,s1*(1-psi),s1*psi,1-s1,
      0,s2*(1-psi),s2*psi,1-s2,
      0,0,s3,1-s3,
      0,0,0,1), nrow=n.states,byrow=TRUE)
  }#t
}#i

#Observation process matrix

S.OBS = array(NA,dim=c(n.states,n.states,totalrel,n.occassions-1))
for(i in 1:totalrel){
  for(t in 1:(n.occassions-1)){
    S.OBS[,,i,t]=matrix(c(
      0,0,0,1,
      0,p,0,1-p,
      0,0,p,1-p,
      0,0,0,1), nrow=n.states,byrow=TRUE)
  }
}

#Simulation code of CH
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
    mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] = rep(1:n.occasions, marked[1:n.occasions,s])
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

# Execute function
sim <- simul.ms(S.STATE, S.OBS, marked)
CH <- sim$CH



#Estimating incidence

n.year=4

#number tested for CWD in occassion t is the total tested from the harvest plus the total tested from the marked individuals

#number tested from harvest assuming all harvest individuals were tested
J.test = N.total*h.rate
J.test2 = N.marked*.15 #assuming the recapture rate of marked indivudals ~15%
J =round(J.test+J.test2)

y.neg = rbinom(n.occassions,J,(1-psi)^n.occassions)













################################################################################
###
### IPM Simulation
###
################################################################################

#From Kery BPA book

phi.juv = .3
phi.f= .51
