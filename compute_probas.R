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

compute_probas_subset <- function(N,n,prev, tau){
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

compute_neff<- function(N, n, prev, tau){
  rho = compute_corr(N, prev, tau)
  return(n/(1+ (n-1)* rho))
}

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

compute_p <- function(N, prev, tau){
  # A = 0
  # for (k in 1:(N-1)){
  #   A = A+ choose(N-1,k) * prev^k*(1-prev)^(N-1-k)* (1-(1-tau)^(k))
  #   #print(A)
  # }
  # return(prev + (1-prev)*A)
  return(1-(1-prev)*(1-tau*prev)^(N-1))
}

compute_corr <- function(N, prev, tau){
  p = compute_p(N, prev, tau)
  
  B = compute_cov(N, prev, tau)
  return(B/(p * (1-p)))
}


compute_cov <- function(N, prev, tau){
  p = compute_p(N, prev, tau)
  A = 0
  B = 0
  for (k in 0:(N-2)){
    A = A+ choose(N-2,k) * prev^k*(1-prev)^(N-2-k)* (1-(1-tau)^(k+1))
    if (k>0){
      B = B+ choose(N-2,k) * prev^k*(1-prev)^(N-2-k)* (1-(1-tau)^k)^2}
  }
  return(prev^2 + 2 * prev * (1-prev) *A+ (1-prev)^2 * B - p^2)
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
