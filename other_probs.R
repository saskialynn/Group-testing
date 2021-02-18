compute_discrepancy <- function(N, prev, tau){
  p1 = compute_probas(N, prev, tau)
  p2 = sapply(0:N, function(x){dbinom(x,N, prev)})
  
}

compute_efficient_sample_size <- function(N, prev, tau, dist="Kolmogorov"){
  p1 = compute_probas(N, prev, tau)
  if (dist=="Hellinger"){
    d = sapply(seq(from=1, to= N, by=1), function(n0){
      p2 = sapply(0:N, function(x){dbinom(x,n0, prev)})
      return(1/sqrt(2)*sqrt(sum((sqrt(p1$p) - sqrt(p2))^2)))
    })
  }else{
    d = sapply(seq(from=1, to= N, by=1), function(n0){
      p2 = sapply(0:N, function(x){dbinom(x,n0, prev)})
      return(max(abs(p1$p - p2)))})
  }
  return(list(N_eff = which.min(d), res= d))
}

# compute_p <- function(N, prev, tau){
#   res = rep(0, N+1)
#   res[1] = (1-prev)^N
#   for (K in 1:N){
#     for (k in 1:K){
#       res[K+1] = res[K+1]  + choose(N,k) * prev^k *(1-prev)^(N-k) * choose(N-k, K-k) * (1-tau)^(k*(N-K)) * (1-(1-tau)^k)^(K-k)
#     }
#   }
#   res_temp = data.frame(k = 0:N,
#                         p = res)
#   return(res_temp)
# }


compute_corr_all_levels2 <- function(N, prevs, taus, p = NULL, B=10000){
  #### Have to be way smarter than this...
  if (is.null(p)){
    p = compute_p_all_levels(N, prevs, taus)
  }
  probs <- factor(apply(sapply(prevs, function(x){rbinom(B,1,x)}),1,function(x)sum(x[1:(N-2)])), levels= 0:(N-2))
  probs = as.data.frame(table(probs))
  probs$Freq = probs$Freq/B
  
  prob = sapply(1:N, function(i){
    sapply((1:N)[-i], function(j){
      probs <- factor(apply(sapply(prevs[-c(i,j)], function(x){rbinom(B,1,x)}),1,sum), levels= 0:(N-2))
      probs = as.data.frame(table(probs))
      probs$Freq = probs$Freq/B
      return(
        sum(apply(sapply(1:B, function(b){
          sapply(1:(N-2), function(K){
            infected = sample((1:N)[-c(i,j)], K, prevs[-c(i,j)], replace=FALSE)
            (1-exp(sum(log(1- taus[infected,i])))) * (1-prevs[i]) * prevs[j] + 
              (1-exp(sum(log(1- taus[infected,j])))) * (1-prevs[j]) * prevs[i]+
              + 
              (1-prevs[i]) * (1-prevs[j]) * (1-exp(sum(log(1- taus[infected,i])))) * (1-exp(sum(log(1- taus[infected,j]))))
          })
        }),1,mean) *probs$Freq[2:(N-1)])  + prevs[i] * prevs[j]
      )})
  })
  return( mean(prob -p^2))
  
  
  #   A = 0
  #   Bb = probs[which(probs$probs ==  0), "Freq"]  * 
  #     mean(sapply(1:10000, function(b){(1-exp(sum(log(1-taus[sample(1:length(taus),1)])))) *  exp(sum(log(1-taus[sample(1:length(taus),N-1)])))}))
  #   for (K in 1:(N-2)){
  #     A = 0
  #     Bb = 0
  #     for(k in 0:K){
  #       if (k>0){
  #         Bb = Bb+ probs[which(probs$probs ==  k), "Freq"]  *
  #           mean(sapply(1:1000, function(b){exp(sum(log(1-taus[sample(1:N^2,(N-2-K) * (k+1))]))) * exp(sum(sapply(1:(K+1-k), function(u){log(1-exp(sum(log(taus[sample(1:N^2, k+1)]))))})))}))
  #         A =  A + probs[which(probs$probs ==  k), "Freq"]   * 
  #           mean(sapply(1:1000, function(b){exp(sum(log(1-taus[sample(1:N^2,(N-2-K) * k)]))) * exp(sum(sapply(1:(K-k+2), function(u){log(1-exp(sum(log(taus[sample(1:N^2, k)]))))})))}))
  #         
  #       }else{
  #         Bb = Bb + probs[which(probs$probs ==  0), "Freq"]  * 
  #           mean(sapply(1:10000, function(b){(1-exp(sum(log(1-taus[sample(1:length(taus),1)])))) *  exp(sum(log(1-taus[sample(1:length(taus),N-2)])))}))
  #       }
  # 
  #     }
  #   }
  #   p0 =  mean(sapply(1:1000, function(x){exp(sum(log(1-prevs[sample(1:N, 2)])))}))
  #   p1 =  mean(sapply(1:1000, function(x){ (1-prevs[sample(1:N, 1)]) *(prevs[sample(1:N, 1)]) }))
  #   p2 =  mean(sapply(1:1000, function(x){exp(sum(log(prevs[sample(1:N, 2)])))}))
  #   return((p0 * A + 2*p1 * B + p2- p^2)/(p*(1-p)))
}



compute_corr_all_levels2 <- function(N, prevs, taus, p = NULL, B=10000){
  #### Have to be way smarter than this...
  if (is.null(p)){
    p = compute_p_all_levels(N, prevs, taus)
  }
  probs <- factor(apply(sapply(prevs, function(x){rbinom(B,1,x)}),1,function(x)sum(x[1:(N-2)])), levels= 0:(N-2))
  probs = as.data.frame(table(probs))
  probs$Freq = probs$Freq/B
  
  prob = sapply(1:N, function(i){
    sapply((1:N)[-i], function(j){
      probs <- factor(apply(sapply(prevs[-c(i,j)], function(x){rbinom(B,1,x)}),1,sum), levels= 0:(N-2))
      probs = as.data.frame(table(probs))
      probs$Freq = probs$Freq/B
      return(
        sum(apply(sapply(1:B, function(b){
          sapply(1:(N-2), function(K){
            infected = sample((1:N)[-c(i,j)], K, prevs[-c(i,j)], replace=FALSE)
            (1-exp(sum(log(1- taus[infected,i])))) * (1-prevs[i]) * prevs[j] + 
              (1-exp(sum(log(1- taus[infected,j])))) * (1-prevs[j]) * prevs[i]+
              + 
              (1-prevs[i]) * (1-prevs[j]) * (1-exp(sum(log(1- taus[infected,i])))) * (1-exp(sum(log(1- taus[infected,j]))))
          })
        }),1,mean) *probs$Freq[2:(N-1)])  + prevs[i] * prevs[j]
      )})
  })
  return( mean(prob -p^2))
}


compute_corr_tau_var <- function(N, prev, tau, B=10000){
  #### Compute the correlation (assumes whole group is tested)
  #### when tau varies
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
  B = 0
  for (k in 0:(N-2)){
    A = A+ choose(N-2,k) * prev^k*(1-prev)^(N-2-k)* (1-exp(sum(log(1-tau[k+1]))))
    if (k>0){
      B = B+ choose(N-2,k) * prev^k*(1-prev)^(N-2-k)* (1-exp(sum(log(1-tau[2:(k+1)]))))^2
    }
  }
  return((prev^2 + 2 * prev * (1-prev) *A+ (1-prev)^2 * B - p^2)/(p*(1-p)))
}


compute_effective_parameters<- function(){
  p1 = compute_probas(n0, prev, tau)
  k = do.call(rbind, lapply(1:n0, function(ne){
    t(sapply(seq(from=0.0001, to=0.10, by=0.001), function(prev_eff){
      p2 = sapply(0:ne, function(x){dbinom(x,ne, prev_eff)})
      p2 = c(p2, rep(0, n0-ne))
      c(prev_eff, ne, max(abs(p1$p - p2)))
    }))}))
  ind = which.min(k[,3])
  return(k[ind,])
}


compute_effective_parameters2 <- function(N, prev, tau, dist="Kolmogorov"){
  p1 = compute_probas(N, prev, tau)
  p_e = compute_p(N, prev,tau)
  rho = compute_corr(N, prev, tau)
  N_e  = N/(1+ (N-1)* rho)
  p2 = c(sapply(0:round(N_e), function(x){dbinom(x,round(N_e), p_e )}), rep(0, N-round(N_e)))
  return(c(p_e, N_e, max(p2-p1$p)))
}



compute_probas_level_all_variations <- function(N, prev, tau, tau_relative_var, prev_relative_var, B=10000){
  
  simulations<- rbindlist(lapply(1:B, function(b){
    tilde_tau = matrix(rnorm(N^2, log(tau/(1-tau)), sd=log(tau_relative_var)/2), nrow = N)
    tilde_tau  = 1/(1+exp(- tilde_tau))
    tilde_prev = rnorm(N, log(prev/(1-prev)), sd=log(prev_relative_var)/2)
    tilde_prev  = 1/(1+exp(- tilde_prev))
    pi_eff = compute_p_all_levels(N,tilde_prev,tilde_tau)
    rho = compute_corr_all_levels(N,tilde_prev,tilde_tau)
    probs <- factor(apply(sapply(tilde_prev, function(x){rbinom(B,1,x)}),1,sum), levels= 0:(N))
    probs = as.data.frame(table(probs))
    probs$Freq = probs$Freq/B
    #res = rep(0, N+1)
    #res[1] = exp(sum(log(1-prevs)))
    
    for (K in 1:N){
      sim = sapply(1:B, function(x){
        infected = sample(1:N, K, prevs, replace = FALSE)
        sum(sapply(1-exp(apply(log(1-as.matrix(tilde_tau[, infected])),1,sum)), function(p){rbinom(1,1,p)}))
      })
      sim =    as.data.frame(table( factor(K + sim, levels= 0:(N))))
      sim$Freq = sim$Freq/B *  probs$Freq[probs$probs == K]
      if (K== 1){
        res_tot = sim 
      }else{
        res_tot$Freq = res_tot$Freq + sim$Freq 
      }
    }
    
    res_tot$Freq[1] = exp(sum(log(1-prevs)))
    res_tot$Freq = res_tot$Freq/sum(res_tot$Freq)
    
    # for (K in 1:(N)){
    #   for(k in 1:(K)){
    #     res[K+1]=   res[K+1] + probs[which(probs$probs ==  k), "Freq"]   * choose(N,k) * choose(N-k,N-k) * exp(sum(log(1-taus[sample(1:N^2,(N-K) * k)]))) * 
    #         exp(sum(sapply(1:(K-k), function(u){log(1-exp(sum(log(1-taus[sample(1:N^2, k)]))))})))
    #     
    #   }
    # }
    res_temp = data.frame(n = 0:N,
                          p = res_tot$Freq,
                          pi_eff = pi_eff,
                          n_eff= N/(1+(N-1) *rho),
                          prob_null_eff = sapply(0:N, function(positives){dbinom(positives, N, pi_eff )}))
    
  }))
  return(simulations %>% group_by(n) %>% summarize_all(mean))
}
