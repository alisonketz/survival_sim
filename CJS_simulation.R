#CJS simulation


# Define parameter values
n.occasions = 5, number of capture occasions
marked = rep(100, n.occasions-1) #annual number of newly marked individuals

phi=rep(0.4,n.occasions-1) #assume survival is .4
p=rep(.25,n.occasions-1) #assume recapture probability is .25


#Define matrices with survival and recapture probabilities
PHI = matrix (phi,ncol=n.occassions, nrow=sum(marked))



