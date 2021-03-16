

compute_pd <- function(N, K, prev){
  sapply(1:K, function(k){choose(N,k) * prev^k * (1-prev)^(N-k) })
}

compute_td <- function(N, K, tau){
  sapply(1:(K), function(k){choose(N-k,K-k) *((1-tau)^k)^(N-K) * (1-(1-tau)^k)^(K-k)})
}

prev_graph_effect <- function(x, B){
  a = rnorm(B, log(x/(1-x)), 1)
  return((1/(1+exp(-a))))
}

prev_subject_effect <- function(a, N){
  a = sapply(x, function(xx){rnorm(N,log(xx/(1-xx)), 1)})
  return((1/(1+exp(-a))))
}

proba_laws <- function(N, prev, tau, 
                       tau_graph_effect=NULL,
                       tau_subject_effect=NULL, 
                       prev_graph_effect=NULL,
                       prev_subject_effect=NULL,
                       B=1000){
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
  if (is.null(prev_graph_effect) == FALSE){
    prevs = prev_graph_effect(prev, B)
    dist_nothing <- sapply(prevs, function(prev){(1-prev)^N})
    dist_prevs = lapply(prevs, function(prev){
      c((1-prev)^N, sapply(1:N, function(K){ compute_pd(N,K, prev)
        }))})  ### list (B length) of K items of length k
  }else{
    ### Let's deal with this latter
        prevs = rep(prev, B)
        dist_nothing <- sapply(prevs, function(prev){(1-prev)^N})
        dist_prevs = lapply(prevs, function(prev){
          c((1-prev)^N, sapply(1:N, function(K){ compute_pd(N,K, prev)
          }))})   ### list (B length) of K items of length k
        #choose(N-k, K-k) * exp( sum(log(1-taus((N-K)*k)))) * exp(sum(sapply(1:(K-k), function(toto){log(1-exp(sum(log(1-taus(k)))))})))
  }
  if (is.null(prev_subject_effect) == FALSE){
    #### the prevalence is small, so we use a Poisson approximation
    dist_prevs = sapply(prevs, function(prev){
      #### generate subject effects
      prev_s = prev_subject_effect(prev, N)
      lambda=  sum(prevs)
      return(sapply(0:N, function(K){exp(-lambda) * lambda^K/factorial(K)}))
    })}

  
  if (is.null(tau_graph_effect) == FALSE){
    taus = tau_graph_effect(tau, B)
    dist_taus = lapply(taus, function(tau){
      sapply(1:N, function(K){ compute_td(N,K, tau)
      })})  ### N x B
  }else{
    taus = rep(tau, B)
    dist_taus = lapply(taus, function(tau){
      sapply(1:N, function(K){ compute_td(N,K, tau) 
      })})
  }
  if (is.null(tau_subject_effect) == FALSE){
    taus_f <- function(b){tau_subject_effect(b, tau=tau)}
    dist_taus = lapply(1:B, function(b){
      sapply(1:N, function(K){
      #### generate subject effects
      sapply(1:(K), function(k){
        taus_v = taus_f((N-K)*K) 
        if (k<K){
          choose(N-k, K-k) * exp( sum(log(1-taus_f((N-K)*k)))) * exp(sum(sapply(1:(K-k),
                                                                              function(j){log(1-exp(sum(log(1-taus_f(k)))))})))  
        }else{
          exp( sum(log(1-taus_f((N-K)*K))))
        }
        
      })
    })})
  }
  
  
  #### Now compute the dot product

    test = sapply(1:B,function(b){
        c(dist_prevs[[b]][[1]],sapply(1:N, function(k){
          sum(dist_prevs[[b]][[(k+1)]] * dist_taus[[b]][[k]])
        }))
      }) 
  for (j in 1:ncol(test)){
    if(sum(test[,j])!=1){
      test[,j] = test[,j]/sum(test[,j])
    }
    
  }
  return(test)
}
  

#### Is this legit?  ## Need to run the simulations






