
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
    
