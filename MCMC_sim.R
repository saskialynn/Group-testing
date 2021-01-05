calc_probs <- function(N, prev, min, max, B=10000){ # N = group size, prev = prevalence
  infector <- numeric(B) # initial infectious individual in network
  infected <- matrix(nrow = B, ncol = N) # secondary infected individuals
  dimnames(infected) <- list(1:B, c(1:N)) # rows = experiments, columns = # positives
  v = prev * rep(1,N)
  proba = 1-exp(sum(log(1 - v))) # probability of somebody being infected
  for (i in 1:B){
    infector[i] <- which(rmultinom(1, size = 1, prob= (v/sum(v))) > 0) # choose initial infectious individual
    if (N >1){ # group size > 1 
      infected[i, infector[i]] <- 1 # infector 
      rest = setdiff(1:N, infector[i]) # susceptible individuals
      infected[i, rest] <- pmin(1, sapply(rest, function(x){rbinom(1,1,runif(1, min, max))}) + sapply(rest, function(x){rbinom(1,1,prev)}))
      # infected[i, rest] <- pmin(1, rbinom(rest, 1, runif(rest, min, max)) + rbinom(rest,1,prev))
      # propagate infection: someone can be infected by the infector individual or from outside the network (prevalence)
      
    }else{
      infected[i] = 1
    }

  }
  
  
  # store probability of sum(x)=k positives for given group size
  temp <- apply(infected, 1, sum) # sum over the rows of infected 
  temp <- factor(temp, levels= 0:N) 
  res_temp = as.data.frame(table(temp))
  res_temp$Freq = res_temp$Freq/B * proba # weight by probability of somebody being infected
  res_temp$Freq[1] = 1 - proba # probability of nobody being infected
  names(res_temp) = c("n", "Freq")
  return(res_temp)
}









