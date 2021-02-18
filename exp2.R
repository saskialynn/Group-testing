a = 0 
N=15

compute_p_all_levels <- function(N, prevs, tau,var_tau=2, B=100000){
  #### when both prev and tau are vectors
  probs <- factor(apply(sapply(prevs, function(x){rbinom(B,1,x)}),1,
                        function(x){sum(x[-sample(1:N,1)])}), levels= 0:(N-1))
  probs = as.data.frame(table(probs))
  probs$Freq = probs$Freq/B
  probs$probs = as.numeric(probs$probs)
  return(
    sum(apply(sapply(1:B, function(b){
      sapply(1:(N-1), function(K){
        (1-exp(sum(log(1- 1/(1+exp(-rnorm(K,log(tau/(1-tau)), log(var_tau)/2)))))))})}),1,mean) *probs$Freq[2:N])* (1-mean(prevs)) + mean(prevs)
  )
  # return(
  #   sum(apply(sapply(1:B, function(b){
  #     sapply(1:(N-1), function(K){
  #       (1-exp(sum(log(1- 1/(1+exp(-rnorm(K,log(tau/(1-tau)), log(var_tau)/2)))))))})}),1,mean) *probs$Freq[2:N])* (1-mean(prevs)) + mean(prevs)
  # )
}


prev  = 0.01
a=0
for (prev in c(0.001, 0.01, 0.1, 0.2)){
for (tau in c( 0, 0.0001, 0.001,0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9) ){
  for (var_tau in seq(from=1, to=6, by=0.5) ){
    print(c(prev,tau))
    start_time <- Sys.time()
    res =rbind(lapply(1:300, function(it)
        {data.frame(p = compute_p_all_levels(N, runif(N, prev*0.5 ,prev*1.5), tau, 
                                             var_tau=var_tau, B=1000),
                      prev = prev,
                      tau=tau,
                      var = var_tau,
                      it = it)}))
      if (a==0){
        res_tot  = res
      }else{
        res_tot = rbind(res_tot, res)
      }
      a  = a+1
      end_time <- Sys.time()
      print(end_time - start_time)
    }
    
  }
}

library(tidyverse)
library(ggplot2)
summary  = res_tot %>%
  group_by(prev,tau, var) %>%
  summarize(mean = mean(p), sd = sd(p), q97=quantile(p,0.975),q2=quantile(p,0.025))

summary$var <- as.factor(summary$var)
ggplot(summary)+
  geom_point(aes(x=tau, y=mean, colour = var)) + 
  geom_ribbon(aes(x=tau, ymin=q2,ymax=q97, fill=var), alpha=0.1 )

