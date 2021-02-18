source("compute_correlations.R")
library(data.table)


compute_probas <- function(N, prev, tau){
  res = rep(0, N+1)
  res[1] = (1-prev)^N
  for (K in 1:N){
    for (k in 1:K){
      res[K+1] = res[K+1]  + choose(N,k) * prev^k *(1-prev)^(N-k) * choose(N-k, K-k) * (1-tau)^(k*(N-K)) * (1-(1-tau)^k)^(K-k)
    }
  }
  res_temp = data.frame(n = 0:N,
                        p = res)
  return(res_temp)
}

compute_probas_tauprev_var <- function(N, prev, tau, alpha=0, alpha_prev=0,
                                            B=10000, mode ="none",
                                       mode_prev="none"){
  #### Compute p with variable tau  and prev (assumes whole group is tested)
  #### when both prev and tau are vectors
  #### INPUT
  #### -------------------------------------------------------------------------
  #### N       : group size
  #### prevs   : prevalences (vector of size N)
  #### tau     : transmissibility within the group 
  #### alpha   : parameter of the variability of transmissibility within the group  (matrix N xN)
  #### mode    : generating mechanism for taus
  #### -------------------------------------------------------------------------
  taus  <- function(k){
    beta = ifelse(alpha > tau, tau,  ifelse(alpha >1-tau, 1-tau, alpha) )
    if(mode == "multiplicative"){
      return(1+ exp(-rnorm(k, log(tau /(1-tau)), sd = log(alpha)/2)))
    }else{
      if(mode == "uniform"){
        return(runif(k, tau - beta, tau + beta ))
      }
      else{
        return(rep(tau,k))
      }
    }
  }
  prevs  <- function(k){
    beta_prev = ifelse(alpha_prev > prev, prev,  ifelse(alpha_prev >1-prev, 1-prev, alpha_prev) )
    if(mode_prev == "multiplicative"){
      return(1+ exp(-rnorm(k, log(prev /(1-prev)), sd = log(alpha_prev)/2)))
    }else{
      if(mode_prev == "uniform"){
        return(runif(k, prev - beta_prev, prev + beta_prev))
      }else{
        return(rep(prev,k))
      }
    }
  }

  
  pi_eff = compute_p_tauprev_var(N, prev, tau, alpha, alpha_prev,
                                 B=10000, mode=mode, mode_prev = mode_prev)
  rho = compute_corr_tauprev_var(N, prev, tau, alpha, alpha_prev, p = pi_eff,
                                 B=10000, mode=mode, mode_prev = mode_prev)
  probs <- factor(sapply(1:B, function(b){ 
    sum(sapply(prevs(N), function(x){rbinom(1,1,x)}))}), levels= 0:(N))
  probs = as.data.frame(table(probs))
  probs$Freq = probs$Freq/B
  
  
  simulations<- rbindlist(lapply(1:B, function(b){
    res = rep(0, N+1)
    res[1] =  probs$Freq[1]
    res[2] = probs$Freq[2] * exp( sum(log(1-taus((N-1)))))
    for (K in 2:N){
      for (k in 1:K){
        res[K+1] = res[K+1]  + 
          probs$Freq[k+1] * choose(N-k, K-k)  * exp( sum(log(1-taus((N-K)*k)))) * exp(sum(sapply(1:(K-k), function(toto){log(1-exp(sum(log(1-taus(k)))))})))
      }
    }
    if (abs(sum(res) -1) > 0.01){
      return("Error computation for the probability")
    }else{
    res = res/sum(res)  #### (numerical imprecision up to order B due to the approximation of B)
    }
    
    res_temp = data.frame(n = 0:N,
                          p = res,
                          rho=rho,
                          pi_eff = pi_eff,
                          n_eff= N/(1+(N-1) *rho),
                          prob_null_eff = sapply(0:N, function(positives){dbinom(positives, N, pi_eff )}))
    
  }))
  return(simulations %>% group_by(n) %>% summarize_all(mean))
}

compute_probas_subset <- function(N,n,prev, tau){
  #### If we only observe part of a group
  res = rep(0, n + 1)
  p = (compute_probas(N, prev, tau))$p
  for (K in 0:n){
    for (k in 0:(N-n)){
      res[K+1] = res[K+1] + p[1 + k + K] *  choose(K+k, K)* choose(N-K-k, n-K)/choose(N, n)
    }
  }
  res_temp = data.frame(n = 0:n,
                        p = res)
  return(res_temp)
}


