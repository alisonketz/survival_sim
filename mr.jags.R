
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



