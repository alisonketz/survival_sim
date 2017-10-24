
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
    
