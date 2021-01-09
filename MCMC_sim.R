calc_probs <- function(N, prev, min, max, B=10000, risk=NULL){ # N = group size, prev = prevalence
  infector <- numeric(B) # initial infectious individual in network
  infected <- matrix(nrow = B, ncol = N) # secondary infected individuals
  dimnames(infected) <- list(1:B, c(1:N)) # rows = experiments, columns = # positives
  v = ifelse(is.null(risk), prev * rep(1,N), v) ### we are assuming everyone has the same probability of becoming infectious
  proba = 1-exp(sum(log(1 - v))) # probability of at least somebody being infected
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


calc_probs_all_homogeneous<- function(N, prev, min, max, B=1e5, risk=NULL){ # N = group size, prev = prevalence
  infector <- numeric(B) # initial infectious individual in network
  infected <- matrix(nrow = B, ncol = N) # secondary infected individuals
  dimnames(infected) <- list(1:B, c(1:N)) # rows = experiments, columns = # positives
  if(is.null(risk)){
    v = prev * rep(1,N)
  }else{
    v = risk
  }### we are assuming everyone has the same probability of becoming infectious
  proba = 1-exp(sum(log(1 - v))) # probability of at least somebody being infected
  res = sapply(1:B, function(b){
      #k = sum(sapply(1:N, function(i){ rbinom(1,1,v[i])} ))
      k = rbinom(1,N, prev)
      add = 0
      #### Step 2: secondary attack
      if(k>0){
        add = sum(sapply(1:N-k, function(j){rbinom(1, 1, 1-exp(sum(log(1-runif(k, min, max)))))}))
      }
    return(k+add)
  } )
  temp <- factor(res, levels= 0:N)
  res_temp = as.data.frame(table(temp))
  res_temp$Freq = res_temp$Freq/B 
  names(res_temp) = c("n", "Freq")
  return(res_temp)
}







