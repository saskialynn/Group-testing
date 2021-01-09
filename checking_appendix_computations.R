#### Checking the approximation we are making


real_proba <- function(G, pi, tau, B=1e6){
  #### Step 1: select the people who are infected
  res = data.frame(K = 0:G,
                   count = rep(0, G+1),
                   approx = rep(0, G+1))
  for (b in 1:B){
    k = rbinom(1, G, pi)
    k0 = rbinom(1, 1, 1-(1-pi)^G)
    add = 0
    add0 = 0
    #### Step 2: secondary attack
    if(k>0){
      add = rbinom(1, G-k, 1-(1-tau)^k)
    }
    res$count[which(res$K == add + k)] = 1/B + res$count[which(res$K == add + k)]
    
    if(k0>0){
      add0 = rbinom(1, G-k0, tau)
    }
    res$approx[which(res$K == add0 + k0)] = 1/B + res$approx[which(res$K == add0 + k0)]
  }
  res$delta= res$count - res$approx
  return(res)
}

compute_delta2 <- function(G, K, pi, tau){
  tilde_tau = 1 - tau
  vec = sapply(0:K, function(k){
    choose(G, k) * pi^k * (1-pi)^(G-k) * choose(G - k, K -k) * (1-tilde_tau^k)^(K-k) * tilde_tau^(k* (G-K))
  })
  if (K>1){
    A = sum(vec[3:(K+1)])
  }else{
    A=  0
  }
  return(A)
}

compute_delta3 <- function(G, K, pi, tau){
  tilde_tau = 1 - tau
  vec = sapply(0:K, function(k){
    choose(G, K) * choose(K, k) *  pi^k *  ((1-pi) *(1-tilde_tau^k))^(K-k) * ((1-pi)  * tilde_tau^k)^ (G-K)
  })
  A = sum(vec)
  return(A)
}

compute_delta <- function(G, K, pi, tau){
  tilde_tau = 1-tau
  vec = sapply(c(1:K), function(k){
    choose(K,k) * pi^k * ( (1-pi)*(1-tilde_tau^k))^(K-k) * ((1-pi) * tilde_tau^k)^(G-K)
  })
  A = choose(G,K) * sum(vec)
  return(A)
}


compute_delta <- function(G, K, pi, alpha){
  tilde_tau = 1 - alpha * pi
  vec = sapply(0:K, function(k){
    choose(G, k) * pi^k * (1-pi)^(G-k) * choose(G - k, K -k) * (1-tilde_tau^k)^(K-k) * tilde_tau^(k* (G-K))
  })
  if (K>1){
    A = sum(vec[3:(K+1)])
  }else{
    A=  0
  }
  return(A)
}


pis = seq(from=0.00001, to= 0.2, by = 0.00001)
alphas = seq(from=1, to= 10, by = 0.01)
df = as.data.frame(
  sapply(pis, FUN=function(pi){
  sapply(alphas, function(alpha){compute_delta2(10, 5, pi, alpha)})
  }))

library(tidyverse)
df%>%ggplot(aes(pi,alpha))+stat_density2d(geom="polygon",
                                          aes(alpha = ..level..),
                                          fill="orangered",color="red4",linetype=2)+
  theme_bw()+scale_x_log10("X-coordinate")+scale_y_continuous("Y-coordinate")

library(reshape2)
test = dcast(df, formula = pi ~ alpha)
heatmap(df, Colv = NA, Rowv = NA, scale="column",  xlab="pi", ylab="alpha", main="heatmap")

