

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
  

pool.size = n
B=1000
sens = sapply(1:N, function(positives){
  dat <- matrix(sample(threshold.ct, positives * B, replace=T), nrow=positives) 
  # sample data uniformly at random 
  # n samples of positives, rearrange into a matrix of positive rows, n columns
  each.conc = -log2(colSums(2^-dat)/pool.size)+ifelse(dilution.vary.index==1,0,
                                                      rnorm(mean=0,sd=1.1,n=ncol(dat)))
  z.index= probit.z.indices[probit.mode.index]
  return(mean(probit[1+(z.index-1)*571+each.conc*10-(lod-35.9)*10,2]))

})


prev_graph_effect <- function(x, B){
  a = rnorm(B, log(x/(1-x)), 3)
  return((1/(1+exp(-a))))
}

tau_graph_effect <- function(x, B){
  return(runif(B, 0.5 *tau, 1.5 * tau ))
}

it =1
for (prev in c(0.001, 0.001, 0.005, 0.01, 0.02, 
               0.03, 0.05, 0.07, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3)){
  for(tau in seq(from=0.00, to=0.7, by=0.005)){
    print(c(prev,tau))
    a = list(proba_laws(N, prev, tau, 
                   tau_graph_effect=NULL,
                   tau_subject_effect=NULL, 
                   prev_graph_effect=NULL,
                   prev_subject_effect=NULL,
                   B=1000),
             proba_laws(N, prev, tau, 
                        tau_graph_effect=tau_graph_effect,
                        tau_subject_effect=NULL, 
                        prev_graph_effect=NULL,
                        prev_subject_effect=NULL,
                        B=1000),
             proba_laws(N, prev, tau, 
                        tau_graph_effect=NULL,
                        tau_subject_effect=NULL, 
                        prev_graph_effect=prev_graph_effect,
                        prev_subject_effect=NULL,
                        B=1000),
             proba_laws(N, prev, tau, 
                        tau_graph_effect=tau_graph_effect,
                        tau_subject_effect=NULL, 
                        prev_graph_effect=prev_graph_effect,
                        prev_subject_effect=NULL,
                        B=1000)
             )
    names = c("Fixed", "Tau Graph Effect","Pi Graph Effect","All Graph Effect")
    for (n in 1:4){
      if (it ==1){
        res = data.frame(x = sapply(1:B, function(b){
          sum(a[[n]][2:(N+1), b] * sens)/sum(a[[n]][2:(N+1), b])}),
          type = names[n])
      }else{
        res = rbind(res,
                    data.frame(x = sapply(1:B, function(b){
          sum(a[[n]][2:(N+1), b] * sens)/sum(a[[n]][2:(N+1), b])}),
          type = names[n]))
      }
      it = it + 1
      
    }
    
  }
  
}


res$
ggplot(res) +
  
  


test = sapply(1:B, function(b){
 sum(a[, b] * c(0,sens))
})

test = sapply(1:B, function(b){
  sum(a[2:(N+1), b] * sens)/sum(a[2:(N+1), b])
})
#### Is this legit?  ## Need to run the simulations






