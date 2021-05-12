# Modified code from Claire's compute_proba_laws.R
##############################
####  Compute proba laws (approximation) with whilst allowing for some variability
####  INPUT
####  --------------------------------------------
####  N, prev, tau         : Size of the group, baseline prevalence and transmissibility
####  tau_graph_effect     : function to sample variables taus (graph effect). NULL if no graph effect.
####  tau_subject_effect   : function to sample variables taus (edge effect). NULL if no edge effect.
####  prev_graph_effect    : function to sample variables prevs (graph effect). NULL if no graph effect.
####  prev_subject_effect  : function to sample variables prevs (subject effect). NULL if no edge effect.
####

# functions to compute probability of K positives in a pool of size N
compute_pd <- function(N, K, prev){ # community infections
  sapply(1:K, function(k){choose(N,k) * prev^k * (1-prev)^(N-k) })
}

compute_td <- function(N, K, tau){ # network transmission
  sapply(1:K, function(k){choose(N-k,K-k) *((1-tau)^k)^(N-K) * (1-(1-tau)^k)^(K-k)})
}
# functions to sample prevalence and tau from prior distributions
prev_graph_effect <- function(x, B){ # how to change this based on the distribution that I fitted?
  a = rnorm(B, log(x/(1-x)), 1)
  return((1/(1+exp(-a))))
}

prev_subject_effect <- function(x, N){ 
  a = sapply(x, function(xx){rnorm(N,log(xx/(1-xx)), 1)})
  return((1/(1+exp(-a))))
}

tau_graph_effect <- function(x, B){ # sample tau from unif[.5*tau, 1.5*tau]
  return(runif(B, 0.5 *x, 1.5 * x ))
}

tau_subject_effect <- function(x, N){ # check this with Claire
  a = sapply(x, function(xx){runif(N, 0.5*xx, 1.5*xx)}) # change to sampling from pareto 
  return(a)
}

######## Probability Laws
proba_laws <- function(N, prev, tau, 
                       tau_graph_effect=NULL,
                       tau_subject_effect=NULL, 
                       prev_graph_effect=NULL,
                       prev_subject_effect=NULL,
                       null_mod=NULL,
                       B=1000){
  
  ## Null Model
  if(is.null(null_mod) == FALSE){
    prevs=rep(prev, B)
    dist_prevs = sapply(prevs, function(prev){
      dbinom(0:N, N, prev)
    })
    return(dist_prevs)
  }else{
    ## Sample prevalence from prior
    if (is.null(prev_graph_effect) == FALSE){ 
      prevs = prev_graph_effect(prev, B) # sample prevalence B times
      #dist_nothing <- sapply(prevs, function(prev){(1-prev)^N}) # nobody infected from community
      dist_prevs = lapply(prevs, function(prev){ # 1:K infected from comm
        c((1-prev)^N, sapply(1:N, function(K){ compute_pd(N,K, prev)
        }))})  ### list (B length) of N+1 items (for 0:N infected), each of length k
      # dist_prevs is the probability of 1:k individuals infected in pool
    }else{ # point value of prev (not sampled)
      prevs = rep(prev, B)
      #dist_nothing <- sapply(prevs, function(prev){(1-prev)^N})
      dist_prevs = lapply(prevs, function(prev){
        c((1-prev)^N, sapply(1:N, function(K){compute_pd(N,K, prev)
        }))})   ### list (B length) of K items of length k
    }
    
    ## REVIEW LATER
    if (is.null(prev_subject_effect) == FALSE){
      #### the prevalence is small, so we use a Poisson approximation to the binomial
      ## large n, small pi -> lambda = n(pi)
      dist_prevs = sapply(prevs, function(prev){
        #### generate subject effects
        prev_s = prev_subject_effect(prev, N) # sample prev
        lambda=  sum(prev_s)
        return(sapply(0:N, function(K){exp(-lambda) * lambda^K/factorial(K)})) #Poisson PMF
      })}
    
    ## Sample tau from prior
    if (is.null(tau_graph_effect) == FALSE){
      taus = tau_graph_effect(tau, B)
      dist_taus = lapply(taus, function(tau){
        sapply(1:N, function(K){ compute_td(N,K, tau) # compute for each tau and value of N
          # MODIFIED: added c((1-tau)^N)
        })})  ### N x B
    }else{
      taus = rep(tau, B)
      dist_taus = lapply(taus, function(tau){
        sapply(1:N, function(K){ compute_td(N,K, tau) 
        })})
    }
    
    ## REVIEW LATER
    if (is.null(tau_subject_effect) == FALSE){
      taus_f <- function(b){tau_subject_effect(x=tau, b)}
      dist_taus = lapply(1:B, function(b){ #for B simulations
        sapply(1:N, function(K){ 
          #### generate subject effects
          sapply(1:(K), function(k){
            # taus_v = taus_f((N-K)*K) 
            if (k<K){
              choose(N-k, K-k) * exp(sum(log(1-taus_f((N-K)*k)))) * exp(sum(sapply(1:(K-k), function(j){log(1-exp(sum(log(1-taus_f(k)))))})))  
            }else{
              exp(sum(log(1-taus_f((N-K)*K))))
            }
            
          })
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