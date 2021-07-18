# Modified code from Claire's compute_proba_laws.R
##############################
####  Compute proba laws (approximation) with whilst allowing for some variability
####  INPUT
####  --------------------------------------------
####  N, prev, tau         : Size of the group, baseline prevalence and transmissibility
####  prev, tau            : mean of beta(prev_alpha, prev_beta) and beta(tau_alpha, tau_beta)
####  tau_graph_effect     : function to sample variables taus (graph effect). NULL if no graph effect.
####  prev_graph_effect    : function to sample variables prevs (graph effect). NULL if no graph effect.
####

# functions to compute probability of K positives in a pool of size N
compute_pd <- function(N, K, prev){ # community infections
  sapply(1:K, function(k){choose(N,k) * prev^k * (1-prev)^(N-k) })
}

compute_td <- function(N, K, tau){ # network transmission
  sapply(1:K, function(k){choose(N-k,K-k) *((1-tau)^k)^(N-K) * (1-(1-tau)^k)^(K-k)})
}

# functions to sample prevalence and tau from prior distributions
prev_graph_effect <- function(alpha, beta, B){ 
  a = rbeta(B, shape1=alpha, shape2=beta)
  return(a)
}

tau_graph_effect <- function(alpha, beta, B){ 
  a = rbeta(B, shape1=alpha, shape2=beta)
  return(a) 
}

######## Probability Laws
proba_laws <- function(N, prev, prev_alpha, prev_beta, 
                       tau, tau_alpha, tau_beta,
                       tau_graph_effect=NULL,
                       prev_graph_effect=NULL,
                       null_mod=NULL,
                       null_mod_prev_graph=NULL,
                       B=1000){
  
  ## Null Model
  if(is.null(null_mod) == FALSE){
    prevs=rep(prev, B)
    dist_prevs = sapply(prevs, function(prev){
      dbinom(0:N, N, prev)
    })
    return(dist_prevs)
  }else if(is.null(null_mod_prev_graph) == FALSE){
    prevs=prev_graph_effect(prev_alpha, prev_beta, B)
    dist_prevs = sapply(prevs, function(prev){
      dbinom(0:N, N, prev)
    })
    return(dist_prevs)
  }else{
    ## Sample prevalence from prior
    if (is.null(prev_graph_effect) == FALSE){ 
      prevs = prev_graph_effect(prev_alpha, prev_beta, B) # sample prevalence B times
      #dist_nothing <- sapply(prevs, function(prev){(1-prev)^N}) # nobody infected from community
      dist_prevs = lapply(prevs, function(prev){ # 1:K infected from comm
        c((1-prev)^N, sapply(1:N, function(K){ compute_pd(N,K, prev)
        }))})  ### list (B length) of N+1 items (for 0:N infected), each of length k
      # dist_prevs is the probability of 1:k individuals infected in pool
    }else{ # point value of prev (not sampled) (mean of beta dist)
      prevs = rep(prev, B)
      dist_prevs = lapply(prevs, function(prev){
        c((1-prev)^N, sapply(1:N, function(K){compute_pd(N,K, prev)
        }))})   ### list (B length) of K items of length k
    }
    
    
    ## Sample tau from prior
    if (is.null(tau_graph_effect) == FALSE){
      taus = tau_graph_effect(tau_alpha, tau_beta, B)
      dist_taus = lapply(taus, function(tau){
        sapply(1:N, function(K){ compute_td(N,K, tau) # compute for each tau and value of N
        })})  ### N x B
    }else{
      taus = rep(tau, B)
      dist_taus = lapply(taus, function(tau){
        sapply(1:N, function(K){ compute_td(N,K, tau) 
        })})
    }
    
    
    #### Now compute the dot product
    test = sapply(1:B,function(b){
      c(dist_prevs[[b]][[1]], sapply(1:N, function(k){
        sum(dist_prevs[[b]][[(k+1)]] * dist_taus[[b]][[k]])
      }))
    }) 
    # make sure sums to 1 (since we are approximating)
    # normalize if needed 
    for (j in 1:ncol(test)){
      if(sum(test[,j])!=1){
        test[,j] = test[,j]/sum(test[,j])
      }
      
    }
    return(test)
  }
}