# nb_group <- 4
# nb_out <- 5
# N <- nb_group + nb_out
# B <- 1000
# Sigma <- matrix(runif(N^2, 0.1, 0.9), nrow = N, ncol = N) # generate random values from .1 to 0.9
# Sigma = 0.5 * (Sigma + t(Sigma)) # make matrix symmetrical 
# diag(Sigma)= rep(1,N) # 1's on diagonal 
# prev <- .02*rep(1:(N))

simulate_infections <- function(nb_group, nb_out, Sigma, prev, B=1000){
  #### 
  ## v is the vector of individual probabilities of becoming infected
  ## Sigma is the matrix with the edge transmission probabilities
  ## P is the remaining number of people to take into account
  ## B is the number of simulations
  ####
  #print(index_group)
  G = nb_group # number of correlated individuals infected by network
  N = nb_group + nb_out 
  P = nb_out # number uncorrelated individuals infected by community
  infector <- numeric(B)

  if (G >=1){    #### if the group is not of size 0
    infected <- matrix(0, nrow = G, ncol = B)
    v = prev[1:G]
    print(1:G)
    print(prev)
    print(v)
    proba = 1-exp(sum(log(1 - prev[1:G]))) #### proba that at least one is the group is infected
    #withProgress(message = paste0('Running ', B, ' simulations for pool size ', N), value = 0, {
    for (i in 1:B){
      print(i)
      # choose starting infected individual (assume that we know that there is at least one infection)
      # transmit infection
      if (G>1){
        infector[i] <- which(rmultinom(1, size = 1, prob= (v/sum(v))) > 0) # choose initial infectious individual
        infected[infector[i], i] <- 1 # infector 
        rest = setdiff(1:G, infector[i]) # susceptible individuals
        infected[rest, i] <- sapply(rest, function(x){rbinom(1,1,Sigma[infector[i], x])}) # infections transmitted in network
      }else{
        infector[i] <- 1
        infected[infector[i], i] <- 1 # infector 
      }
      
    }
    #})
    temp <- apply(infected, 2, sum) # column sum
    temp <- factor(temp, levels= 0:N)
    res_temp = as.data.frame(table(temp))
    res_temp$Freq = res_temp$Freq/B * proba # weight by probability of somebody being infected
    res_temp$Freq[1] = 1 - proba # probability of nobody being infected
    names(res_temp) = c("n", "Freq")
  }else{
    res_temp = data.frame(n = 0:N, Freq=c(1, rep(0,N)))
  }
  
  if (P >=1){
    temp2 <- sapply(1:B,function(x){sum(sapply(prev[(G+1):N], function(x) rbinom(1,1,x)))}) # total number infections transmitted from community
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
   
   for (l in 0:N){ # loop over all people in group (N total)
     for (k in 0:l){ 
        res[l + 1, 2] = res[l + 1, 2] +  as.numeric(res_temp[which(res_temp$n == k), "Freq"]) * as.numeric(res_temp2[which(res_temp2$n == l- k), "Freq"])
      }
    }
  return(res)
}
