

compute_pd <- function(N, K, prev){
  sapply(1:K, function(k){choose(N,k) * prev^k * (1-prev)^(N-k) })
}

compute_td <- function(N, K, tau){
  sapply(1:K, function(k){choose(N-k,K-k) *(1-tau)^(k*(N-K)) * (1-(1-tau)^k)^(K-k)})

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
  
  if (is.null(prev_graph_effect) == FALSE){
    prevs = prev_graph_effect(prev, B)
    dist_nothing <- sapply(prevs, function(prev){(1-prev)^N})
    dist_prevs = lapply(prevs, function(prev){
      sapply(1:N, function(K){ compute_pd(N,K, prev) 
        })})  ### N x B
  }else{
    prevs = rep(prev, B)
    dist_nothing <- sapply(prevs, function(prev){(1-prev)^N})
    dist_prevs = matrix(sapply(1:N, function(K){ compute_pd(N,K, prev)
    }), ncol=1)
    
    #choose(N-k, K-k) * exp( sum(log(1-taus((N-K)*k)))) * exp(sum(sapply(1:(K-k), function(toto){log(1-exp(sum(log(1-taus(k)))))})))
  }

  
  ### Let's deal with this latter
  if (is.null(prev_subject_effect) == FALSE){
    #### the prevalence is small, so we use a Poisson approximation
    dist_prevs = sapply(prevs, function(prev){
      #### generate subject effects
      prev_s = prev_subject_effect(prev, N)
      lambda=  sum(prevs)
      return(sapply(1:N, function(K){exp(-lambda) * lambda^K/factorial(K)}))
  })}

  
  if (is.null(tau_graph_effect) == FALSE){
    taus = tau_graph_effect(tau, B)
    dist_taus = lapply(taus, function(tau){
      sapply(1:N, function(K){ compute_td(N,K, tau)
      })})  ### N x B
  }else{
    taus = rep(tau, B)
    dist_taus = matrix(sapply(1:N, function(K){ compute_td(N,K, tau)
    }), ncol=1)
  }
  if (is.null(tau_subject_effect) == FALSE){
    dist_prevs = sapply(taus, function(prev){
      #### generate subject effects
      tau_s = tau_subject_effect(tau, N)
      return(sapply(1:N, function(K){exp(-lambda) * lambda^K/factorial(K)}))
    })}
  
  test = sapply(1:B,function(b){
    sapply(1:N, function(k){
      sum(dist_prevs[[b]][[k]] * dist_taus[[b]][[k]])
    })
  }) 
}


sens = sapply(1:N, function(positives){
  dat <- matrix(sample(threshold.ct, positives * n, replace=T), nrow=positives) 
  # sample data uniformly at random 
  # n samples of positives, rearrange into a matrix of positive rows, n columns
  each.conc = -log2(colSums(2^-dat)/pool.size)+ifelse(dilution.vary.index==1,0,
                                                      rnorm(mean=0,sd=1.1,n=ncol(dat)))
  z.index= probit.z.indices[probit.mode.index]
  mean(probit[1+(z.index-1)*571+each.conc*10-(lod-35.9)*10,2])

})

#### Is this legit?  ## Need to run the simulations

df = data.frame(sens = t(sens)%*% test/apply(test,2,sum))



(1-prev)^N + sum(sapply(1:N, function(K){sum(compute_pd(N,K, prev) * (compute_td(N,K, tau)))}))
