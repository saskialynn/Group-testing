calc_probs <- function(N, prev, min, max, B=1000){ # N = group size, prev = prevalence
  infector <- numeric(B) # initial infectious individual in network
  infected <- matrix(nrow = B, ncol = N) # secondary infected individuals
  dimnames(infected) <- list(1:B, c(1:N))
  v = prev * rep(1,N)
  proba = 1-exp(sum(log(1 - v)))
  for (i in 1:B){
    # generate matrix of transmission probabilities 
    # household transmission- based on literature: 
    # Secondary infection rate = 53% (95% confidence interval [CI] = 46%–60%)
    #Sigma <- matrix(runif(N^2, min, max), nrow = N, ncol = N) # generate random values from .1 to 0.9
    infector[i] <- which(rmultinom(1, size = 1, prob= (v/sum(v))) > 0) 
    if (N >1){
      infected[i, infector[i]] <- 1 # infector 
      rest = setdiff(1:N, infector[i])
      infected[i, rest] <- sapply(rest, function(x){rbinom(1,1,runif(N-1, min, max))})
    }else{
      infected[i] = 1
    }

  }
  
  
  # store probability of sum(x)=k positives for given group size
  temp <- apply(infected, 1, sum)
  temp <- factor(temp, levels= 0:N)
  res_temp = as.data.frame(table(temp))
  res_temp$Freq = res_temp$Freq/B * proba
  res_temp$Freq[1] = 1 - proba
  names(res_temp) = c("n", "Freq")
  return(res_temp)
}