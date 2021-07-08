# nb_group <- 4
# nb_out <- 3
# N <- nb_group + nb_out
# B <- 1000
# # Sigma <- matrix(runif(N^2, 0.1, 0.2), nrow = N, ncol = N) # generate random values from .1 to 0.9
# # Sigma = 0.5 * (Sigma + t(Sigma)) # make matrix symmetrical
# # diag(Sigma)= rep(1,N) # 1's on diagonal
# prev <- runif(N, 0, (2*.03))
# min <- 0.1
# max <- 0.2

simulate_infections <- function(nb_group, nb_out, min,max, prev, B=1000){
  #### 
  ## prev = vector of individual probabilities of becoming infected, length N
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
      if(G==1){ # single "correlated" individual
        k = rbinom(1, 1, prev[1]) # bernoulli trial, prob = prev, to decide if individual is infected 
      }
      if (G>1){ # infections in correlated network
        Sigma <- matrix(runif(G^2, min, max), nrow = G, ncol = G) # generate random values for transmission probs
        Sigma = 0.5 * (Sigma + t(Sigma)) # make matrix symmetrical
        diag(Sigma)= rep(1,G) # 1's on diagonal
        
        # tau <- 0.5 # point estimate of the transmission probability
        # r <- (2*tau) -1 # transform tau to between (-1,1)
        # z <- .5*(log(1+r)-log(1-r)) # Fisher Z-transformation
        # G <- 5
        # r_norm <- rnorm(G, mean = z, sd = 1) # normal prior
        # trans_prob <- (exp(2 * r_norm) - 1) / (exp(2 * r_norm) + 1)
        # trans_prob1 <- .5*(trans_prob+1)
        
        
        #### Step 1: choose initital infectious individuals
        infected = sapply(1:G, function(i){rbinom(1, 1, prev[i])}) # for individual i, decide if infected based on prev[i]
        k = sum(infected) # total number initial infected individuals in network
        #print(infected)
        #### Step 2: secondary attack
        if((k>0) &  (k<G)){
            #add = sum(sapply(1:N-k, function(j){rbinom(1, 1, 1-exp(sum(log(1-runif(k, min, max)))))}))
            add = sum(sapply(which(infected==0), function(j){rbinom(1, 1, 1-exp(sum(log(1-Sigma[j,which(infected==1)]))))}))
            # transmit infection from initial infectious individuals to susceptible individuals
            # p(infected by at least one person) = 1 - p(not being infected by any infectious people)
            # count total number of secondary infections (sum())
        }
      }
      if (P>0){ # infections in uncorrelated community
          add2 = sum(sapply((G+1):N, function(i){rbinom(1, 1, prev[i])}))
      }
      return(k + add + add2) # total number of infections in that simulation
  } )
  temp <- factor(res, levels= 0:N) 
  res_temp = as.data.frame(table(temp)) # count frequency of k positives
  res_temp$Freq = res_temp$Freq/B 
  names(res_temp) = c("n", "Freq")
  return(res_temp)
}










