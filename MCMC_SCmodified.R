# nb_group <- 4
# nb_out <- 3
# N <- nb_group + nb_out
# B <- 1000
# Sigma <- matrix(runif(N^2, 0.1, 0.2), nrow = N, ncol = N) # generate random values from .1 to 0.9
# Sigma = 0.5 * (Sigma + t(Sigma)) # make matrix symmetrical
# diag(Sigma)= rep(1,N) # 1's on diagonal
# prev <- runif(N, 0, (2*.03))


simulate_infections <- function(nb_group, nb_out, min,max, prev, B=1000){
  #### 
  ## v = vector of individual probabilities of becoming infected
  ## Sigma = matrix with the edge transmission probabilities
  ## nb_group = number of correlated individuals in network (nb_group = G)
  ## nb_out = remaining number of people to take into account (uncorrelated individuals) (nb_out = P)
  ## B = number of simulations
  ####
  #print(index_group)
  G = nb_group # number of correlated individuals infected by network
  N = nb_group + nb_out 
  P = nb_out # number uncorrelated individuals infected by community
  
  res = sapply(1:B, function(b){
      #k = sum(sapply(1:N, function(i){ rbinom(1,1,v[i])} ))

      add = 0
      add2 = 0
      k=0
      if(G==1){
        k = rbinom(1, 1, prev[1])
      }
      if (G>1){
        Sigma <- matrix(runif(G^2, min, max), nrow = G, ncol = G) # generate random values for transmission probs
        Sigma = 0.5 * (Sigma + t(Sigma)) # make matrix symmetrical
        diag(Sigma)= rep(1,G) # 1's on diagonal
        infected = sapply(1:G, function(i){rbinom(1, 1, prev[i])})
        k = sum(infected)
        #print(infected)
        #### Step 2: secondary attack
        if((k>0) &  (k<G)){
            #add = sum(sapply(1:N-k, function(j){rbinom(1, 1, 1-exp(sum(log(1-runif(k, min, max)))))}))
            add = sum(sapply(which(infected==0), function(j){rbinom(1, 1, 1-exp(sum(log(1-Sigma[j,which(infected==1)]))))}))
        }
      }
      if (P>0){
          add2 = sum(sapply((G+1):N, function(i){rbinom(1, 1, prev[i])}))
      }
      return(k + add + add2)
  } )
  temp <- factor(res, levels= 0:N)
  res_temp = as.data.frame(table(temp))
  res_temp$Freq = res_temp$Freq/B
  names(res_temp) = c("n", "Freq")
  return(res_temp)
}










