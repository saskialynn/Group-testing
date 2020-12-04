simulate_infections <- function(v, Sigma, prev, B=1000){
  #### 
  ## v is the vector of individual probabilities of becoming infected
  ## Sigma is the matrix with the edge transmission probabilities
  ## P is the remaining number of people to take into account
  ## B is the number of simulations
  ####
  
  G = length(v)
  P = length(prev)
  infector <- numeric(B)
  infected <- matrix(nrow = G, ncol = B)
  if (G >=1){
    for (i in 1:B){
      # choose starting infected individual (assume that we know that there is at least one infection)
      # transmit infection
      infector[i] <- which(rmultinom(1, size = 1, prob= (v/sum(v))) > 0) 
      infected[infector[i], i] <- 1 # infector 
      rest = setdiff(1:G, infector[i])
      infected[-infector[i], i] <- sapply(rest, function(x){rbinom(1,1,Sigma[infector[i], x])})
      #infected[i,-infector[i]] <- unlist(lapply(infect_contact$transmission_prob, function(x) x > runif(1, min = 0, max = 1)))
      # infected if contact rate > uniform(0,1)
      proba = 1-exp(sum(log(1-v*prev_tot)))
    }
    if (P>0){
      #### This is not correct
      temp <- apply(infected, 2, sum) + sapply(1:B,function(x){sum(sapply(prev, function(x) rbinom(1,1,x)))})
    }else{
      temp <- apply(infected, 2, sum)
    }
    
    temp <- factor(temp, levels= 1:(G+P))
  }else{
    temp <- sapply(1:B,function(x){sum(sapply(prev, function(x) rbinom(1,1,x)))})
    temp <- factor(temp, levels= 1:(G+P))
  }

  return(as.data.frame(table(temp)/B))
}