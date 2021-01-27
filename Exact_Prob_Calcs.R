# nb_group <- 1
# nb_out <- 1
# N <- nb_group + nb_out
# B <- 1000
# Sigma <- matrix(runif(N^2, 0.1, 0.9), nrow = N, ncol = N) # generate random values from .1 to 0.9
# Sigma = 0.5 * (Sigma + t(Sigma)) # make matrix symmetrical
# diag(Sigma)= rep(1,N) # 1's on diagonal
# prev <- .02*rep(1:(N))

exact_probs <- function(group_size, tau, pop_prev){
  #####
  # group_size = total number of people in group
  # pop_prev = transmission prob in community (prevalence)
  # tau = transmission prob in network
  # group_size <- 10
  # tau <- 0.2
  # pop_prev <- 0.01
  
  N <- group_size
  res <- data.frame(n = 0:N,
                    Freq = rep(0, N+1))


  # probability of L positives in a group of size N
  res[1, 2] <- (1-pop_prev)^N
  for (L in 1:N){ # loop over total positives possible (1:group size)
    for (k in 1:L){ # loop over number of people infected outside
        prob_infect_out <- choose(N,k)*(pop_prev^k)*(1-pop_prev)^(N-k) # k people infected outside network 
        prob_infect_in <- choose((N-k), (L-k))*((1-(1-tau)^k)^(L-k))*((1-tau)^k)^(N-L) # L-k people infected in network 
        res[L + 1, 2] = res[L + 1, 2] + (prob_infect_out*prob_infect_in)
        # res[L + 1, 2] = res[L + 1, 2] + ((choose(N,k)*(pop_prev^k)*(1-pop_prev)^(N-k))*(choose((N-k), (L-k))*((1-(1-tau)^k)^(L-k))*((1-tau)^k)^(N-k-L)))
      }
  }
  return(res)
}

# test <- exact_probs(1, .2, .1)
# sum(test$Freq)



