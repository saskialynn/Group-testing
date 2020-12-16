simulate_infections <- function(nb_group, nb_out, Sigma, prev, B=1000){
  #### 
  ## v is the vector of individual probabilities of becoming infected
  ## Sigma is the matrix with the edge transmission probabilities
  ## P is the remaining number of people to take into account
  ## B is the number of simulations
  ####
  #print(index_group)
  G = nb_group
  N = nb_group + nb_out 
  P = nb_out
  infector <- numeric(B)

  if (G >=1){    #### if the group is not of size 0
    infected <- matrix(nrow = G, ncol = B)
    v = prev[1:G]
    print(1:G)
    print(prev)
    print(v)
    proba = 1-exp(sum(log(1 - prev[1:G]))) #### proba that one is the group is infected
    #withProgress(message = paste0('Running ', B, ' simulations for pool size ', N), value = 0, {
    for (i in 1:B){
      # choose starting infected individual (assume that we know that there is at least one infection)
      # transmit infection
      infector[i] <- which(rmultinom(1, size = 1, prob= (v/sum(v))) > 0) 
      infected[infector[i], i] <- 1 # infector 
      rest = setdiff(1:G, infector[i])
      infected[rest, i] <- sapply(rest, function(x){rbinom(1,1,Sigma[infector[i], x])})
      #incProgress(1/B, detail = paste("*"))
    }
    #})
    temp <- apply(infected, 2, sum)
    temp <- factor(temp, levels= 0:N)
    res_temp = as.data.frame(table(temp))
    res_temp$Freq = res_temp$Freq/B * proba
    res_temp$Freq[1] = 1 - proba
    names(res_temp) = c("n", "Freq")
  }else{
    res_temp = data.frame(n = 0:N, Freq=c(1, rep(0,N)))
  }
  
  if (P >=1){
    temp2 <- sapply(1:B,function(x){sum(sapply(prev[G+1:N], function(x) rbinom(1,1,x)))})
    temp2 <- factor(temp2, levels= 0:(N))
    res_temp2 = as.data.frame(table(temp2))
    res_temp2$Freq = res_temp2$Freq/B
    names(res_temp2) = c("n", "Freq")
    print(res_temp2)
  }else{
    res_temp2 = data.frame(n = 0:N,
                           Freq=c(1, rep(0,N)))
  }
    #### Step 3: Now add the rest of the people that do not necessarily belong to the group
   res <- data.frame(n = 0:N,
                      Freq = rep(0, N+1))
   #print(res_temp)
   #print(res_temp2)
   
   for (l in 0:N){
     for (k in 0:l){
        res[l + 1, 2] = res[l + 1, 2] +  as.numeric(res_temp[which(res_temp$n == k), "Freq"]) * as.numeric(res_temp2[which(res_temp2$n == l- k), "Freq"])
      }
    }
  return(res)
}
