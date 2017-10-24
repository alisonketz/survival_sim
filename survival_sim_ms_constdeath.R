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


############################################################################################################################################################
###
### Multistate Mark Recapture with single disease state
###
############################################################################################################################################################
# 
# s1= .6
# s2 = .3
# psi = .2 #transition probability from susceptible (neg) to infected (pos)
# p1 = .25
# p2 = .25
# n.occasions = 5
# n.states = 3
# n.obs = 3
# marked = matrix(NA,ncol=n.states,nrow = n.occasions)
# marked[,1] = rep(100,n.occasions) #marked negative
# marked[,2] = rep(30,n.occasions) #marked positive
# marked[,3] = rep(0,n.occasions) #marked dead
# 
# 
# # define matrices with survival, transition, and recapture probabilites
# # this is a 4 dimensional matrix
# # dimension 1: state of 'departure'
# # dimension 2: state of 'arrival'
# # dimension 3: individual
# # dimension 4: time
# 
# #state process matrix
# totrel = sum(marked)*(n.occasions-1)
# PSI.STATE = array(NA,dim=c(n.states,n.states,totrel,n.occasions-1))
# for(i in 1:totrel){
#   for(t in 1:(n.occasions-1)){
#     PSI.STATE[,,i,t] = matrix(c(s1*(1-psi),s1*psi,1-s1,0,s2,1-s2,0,0,1),nrow=n.states,byrow=TRUE)
#   }
# }
# 
# p1=c(.25,0,0,0,1)
# 
# 
# #observation process matrix
# PSI.OBS = array(NA,dim=c(n.states,n.obs,totrel,n.occasions-1))
# for(i in 1:totrel){
#   for(t in 1:(n.occasions-1)){
#     PSI.OBS[,,i,t] = matrix(c(p1,0,1-p1,0,p2,1-p2,0,0,1),nrow=n.states,byrow=TRUE)
#   }
# }
# 
# # Define function to simulate multistate capture-recapture data
# simul.ms <- function(PSI.STATE, PSI.OBS, marked, unobservable = NA){
#   # Unobservable: number of state that is unobservable
#   n.occasions <- dim(PSI.STATE)[4] + 1
#   CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
#   # Define a vector with the occasion of marking
#   mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked))
#   g <- colSums(marked)
#   for (s in 1:dim(PSI.STATE)[1]){
#     if (g[s]==0) next  # To avoid error message if nothing to replace
#     mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <-
#       rep(1:n.occasions, marked[1:n.occasions,s])
#   } #s
#   for (i in 1:sum(marked)){
#     for (s in 1:dim(PSI.STATE)[1]){
#       if (mark.occ[i,s]==0) next
#       first <- mark.occ[i,s]
#       CH[i,first] <- s
#       CH.TRUE[i,first] <- s
#     } #s
#     for (t in (first+1):n.occasions){
#       # Multinomial trials for state transitions
#       if (first==n.occasions) next
#       state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
#       CH.TRUE[i,t] <- state
#       # Multinomial trials for observation process
#       event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
#       CH[i,t] <- event
#     } #t
#   } #i
#   # Replace the NA and the highest state number (dead) in the file by 0
#   CH[is.na(CH)] <- 0
#   CH[CH==dim(PSI.STATE)[1]] <- 0
#   CH[CH==unobservable] <- 0
#   id <- numeric(0)
#   for (i in 1:dim(CH)[1]){
#     z <- min(which(CH[i,]!=0))
#     ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id))
#   }
#   return(list(CH=CH[-id,], CH.TRUE=CH.TRUE[-id,]))
#   # CH: capture histories to be used
#   # CH.TRUE: capture histories with perfect observation
# }
# 
# #execute
# 
# ms.sim =  simul.ms(PSI.STATE,PSI.OBS,marked)
# CH = ms.sim$CH
# 
# #compute vector with occasion of first capture
# get.first=function(x)min(which(x!=0))
# f = apply(CH,1,get.first)
# 
# # Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# # 1 = seen alive and neg, 2 = seen alive and pos, 3 = not seen
# rCH <- CH          # Recoded CH
# rCH[rCH==0] <- 3
# 
# # Specify model in JAGS language
# sink("ms.disease.jags.R")
# cat("
# model {
# 
# # -------------------------------------------------
# # Parameters
# # s1: survival probability neg
# # s2: survival probability pos
# # psi: transition probability from neg to pos
# # p1: recapture probability for neg
# # p2: recapture probability for pos
# # -------------------------------------------------
# 
# # -------------------------------------------------
# # States (S):
# # 1 alive at A
# # 2 alive at B
# # 3 dead
# # Observations (O):  
# # 1 seen neg 
# # 2 seen pos
# # 3 not seen
# # -------------------------------------------------
# 
#   # Priors and constraints
#   for (t in 1:(n.occasions-1)){
#      s1[t] <- mean.s[1]
#      s2[t] <- mean.s[2]
#      psi[t] <- mean.psi
#      p1[t] <- mean.p[1]
#      p2[t] <- mean.p[2]
#    }
#   for (u in 1:2){
#      mean.s[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
#      mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
#   }
#   mean.psi ~ dunif(0, 1)    # Priors for mean transitions
# 
#   # Define state-transition and observation matrices
#    for (i in 1:nind){  
#      # Define probabilities of state S(t+1) given S(t)
#      for (t in f[i]:(n.occasions-1)){
#       ps[1,i,t,1] <- s1[t] * (1-psi[t])
#       ps[1,i,t,2] <- s1[t] * psi[t]
#       ps[1,i,t,3] <- 1-s1[t]
#       ps[2,i,t,1] <- 0 
#       ps[2,i,t,2] <- s2[t]
#       ps[2,i,t,3] <- 1-s2[t]
#       ps[3,i,t,1] <- 0
#       ps[3,i,t,2] <- 0
#       ps[3,i,t,3] <- 1
#   # Define probabilities of O(t) given S(t)
#       po[1,i,t,1] <- p1[t]
#       po[1,i,t,2] <- 0
#       po[1,i,t,3] <- 1-p1[t]
#       po[2,i,t,1] <- 0
#       po[2,i,t,2] <- p2[t]
#       po[2,i,t,3] <- 1-p2[t]
#       po[3,i,t,1] <- 0
#       po[3,i,t,2] <- 0
#       po[3,i,t,3] <- 1
#       } #t
#    } #i
#   # Likelihood 
#   for (i in 1:nind){
#      # Define latent state at first capture
#      z[i,f[i]] <- y[i,f[i]]
#      for (t in (f[i]+1):n.occasions){
#         # State process: draw S(t) given S(t-1)
#         z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
#         # Observation process: draw O(t) given S(t)
#         y[i,t] ~ dcat(po[z[i,t], i, t-1,])
#       } #t
#    } #i
# }
# ",fill = TRUE)
# sink()
# 
# # Function to create known latent states z
# known.state.ms <- function(ms, notseen){
#   # notseen: label for ?not seen?
#   state <- ms
#   state[state==notseen] <- NA
#   for (i in 1:dim(ms)[1]){
#     m <- min(which(!is.na(state[i,])))
#     state[i,m] <- NA
#   }
#   return(state)
# }
# 
# z.known = known.state.ms(rCH, 3)
# 
# # Function to create initial values for unknown z
# ms.init.z <- function(ch, f, z.known){
#   for(i in 1:dim(ch)[1]){
#     for(t in 1:n.occasions){
#       if(ch[i,t]==3)ch[i,t]=NA
#     }
#   }
#   
#   for(i in 1:dim(ch)[1]){
#     for(t in f[i]:n.occasions){
#       if(is.na(ch[i,t]))ch[i,t]=ch[i,t-1]
#     }
#   } 
#   for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
#   for(i in 1:dim(ch)[1]){
#     for(t in 1:n.occasions){
#       if(!is.na(z.known[i,t]))ch[i,t]=NA
#     }
#   }  
#   return(ch)
# }
# 
# 
# z.init = ms.init.z(rCH, f,z.known)
# 
# # Bundle data
# jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = z.known)
# 
# # Initial values
# inits <- function(){list(mean.s = runif(2, 0, 1), mean.psi = runif(1, 0, 1), mean.p = runif(2, 0, 1), z = z.init)}  
# 
# # Parameters monitored
# parameters <- c("mean.s", "mean.psi", "mean.p")
# 
# # MCMC settings
# ni = 10000
# nt = 1
# nb = 5000
# nc = 3
# 
# # Call JAGS from R (BRT 8 min)
# ms.cwd <- jags(jags.data, inits, parameters, "ms.disease.jags.R", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
# print(ms.cwd, digits = 3)
# 
# traceplot(ms.cwd)
# 
# ms.cwd.fit=as.mcmc(ms.cwd)
# gelman.diag(ms.cwd.fit)
# 
# dim(ms.cwd.fit[[1]])
# 
# print(ms.cwd)
# 

