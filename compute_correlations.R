compute_neff<- function(N, n, prev, tau){
  #### Compute the effective sample size
  #### INPUT
  #### -------------------------------------------------------------------------
  #### N       : network size
  #### n       : pool size
  #### prev    : prevalence
  #### tau     : transmissibility within the group 
  #### -------------------------------------------------------------------------
  rho = compute_corr(N, prev, tau)
  return(n/(1 + (n-1)* rho))
}



compute_p_tau_var <- function(N, prev, tau, alpha, mode="multiplicative"){
  #### Compute the effective sample size with variable tau (assumes whole group is tested)
  #### INPUT
  #### -------------------------------------------------------------------------
  #### N       : group size
  #### prev    : prevalence
  #### tau     : transmissibility within the group 
  #### alpha   : parameter of the variability of transmissibility within the group  
  #### mode    : generating mechanism for taus
  #### -------------------------------------------------------------------------
  if(mode == "multiplicative"){
    taus =1/(1+ exp(-rnorm(N, log(tau/(1-tau)), sd = log(alpha)/2)))
  }else{
    beta = ifelse(alpha>tau, tau,  ifelse(alpha>1-tau, 1-tau, alpha) )
    taus = runif(N, tau-beta, tau+beta)
  }
  A = 0
  for (k in 1:(N-1)){
    A = A+ choose(N-1,k) * prev^k*(1-prev)^(N-1-k)* (1-exp(sum(log(1-tau[1:k]))))
  }
  return(prev + (1-prev)*A)
}

compute_p_tauprev_var <- function(N, prev, tau, alpha=0, alpha_prev=0, B=10000, 
                                  mode="none",mode_prev="none"){
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
      return(1/(1+ exp(-rnorm(k, log(tau /(1-tau)), sd = log(alpha)/2))))
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
      return(1/(1+ exp(-rnorm(k, log(prev /(1-prev)), sd = log(alpha_prev)/2))))
    }else{
      if(mode_prev == "uniform"){
        return(runif(k, prev - beta_prev, prev + beta_prev))
      }else{
        return(rep(prev,k))
      }
    }
  }
  
  
  probs <- factor(sapply(1:B, function(b){ 
    sum(sapply(prevs(N-1), function(x){rbinom(1,1,x)}))}), levels= 0:(N-1))
  probs = as.data.frame(table(probs))
  probs$Freq = probs$Freq/B
    
  return(
      sum(apply(sapply(1:B, function(b){
        sapply(1:(N-1), function(K){
          (1-exp(sum(log(1- taus(K)))))})}),1,mean) *probs$Freq[2:N])* (1-mean(prevs(B))) + mean(prevs(B))
  )
  
}


compute_corr <- function(N, prev, tau){
  #### Compute the correlation (assumes whole group is tested)
  #### when both prev and tau are vectors
  #### INPUT
  #### -------------------------------------------------------------------------
  #### N       : group size
  #### prev   : prevalence (vector of size N)
  #### tau     : transmissibility within the group 
  p = compute_p(N, prev, tau)
  B = compute_joint(N, prev, tau)
  return((B-p^2)/(p * (1-p)))
}


compute_joint <- function(N, prev, tau){
  #### Compute the joint probability of two variables (assumes whole group is tested)
  #### when both prev and tau are vectors
  #### INPUT
  #### -------------------------------------------------------------------------
  #### N       : group size
  #### prev   : prevalence (vector of size N)
  #### tau     : transmissibility within the group 
  A = 0
  B = 0
  for (k in 0:(N-2)){
    A = A+ choose(N-2,k) * prev^k*(1-prev)^(N-2-k)* (1-(1-tau)^(k+1))
    if (k>0){
      B = B+ choose(N-2,k) * prev^k*(1-prev)^(N-2-k)* (1-(1-tau)^k)^2}
  }
  return(prev^2 + 2 * prev * (1-prev) *A+ (1-prev)^2 * B)
}





compute_corr_tauprev_var <- function(N, prev, tau, alpha, alpha_prev, p = NULL,
                                     B=10000, mode="multiplicative",mode_prev="none"){
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
      return(1/(1+ exp(-rnorm(k, log(tau /(1-tau)), sd = log(alpha)/2))))
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
      return(1/(1+ exp(-rnorm(k, log(prev /(1-prev)), sd = log(alpha_prev)/2))))
      }else{
        if(mode_prev == "uniform"){
          return(runif(k, prev - beta_prev, prev + beta_prev))
        }else{
          return(rep(prev,k))
        }
      }
  }
  
  if (is.null(p)){
    p = compute_p_tauprev_var(N, prev, tau, alpha=alpha, alpha_prev = alpha_prev,
                              B=B, mode=mode, mode_prev = mode_prev)
  }
  
  probs <- factor(sapply(1:B, function(b){ 
    sum(sapply(prevs(N-2), function(x){rbinom(1,1,x)}))}), levels= 0:(N-2))
  probs = as.data.frame(table(probs))
  probs$Freq = probs$Freq/B
  
  pi = mean(prevs(B))
  A = sum(apply(sapply(1:B, function(b){
    log_taus = log(1- taus(2*N))
    a = sapply(1:(N-2), function(K){
      (1-exp(sum(log_taus[1:(K+1)])))})
    b = sapply(1:(N-2), function(K){
      (1-exp(sum(log_taus[1:(K)])))*(1-exp(sum(log_taus[(N+1):(N+K)])))})
    
    return(2* a * (1-pi)  *pi  + b *(1- pi)^2)
  }), 1, mean) * probs$Freq[2:(N-1)] )
  return( (A  + pi^2 - p^2)/(p*(1-p)))
}


