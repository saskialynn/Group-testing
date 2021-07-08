

a = 0 
N=15
for (prev in c(0.0001, 0.001, 0.005, 0.1, 0.2, 0.3, 0.1, 0.15, 0.2)){
  for (tau in c(0, 0.0001, 0.001,0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9) ){
    print(c(prev,tau))
    for (it in 1:500){
      res =data.frame(p = compute_p_all_levels(N, runif(N, prev*0.5 ,prev*1.5), tau, alpha = 0.1, var_tau=2, B=10000),
                      prev = prev,
                      tau=tau,
                      it = it)
      if (a==0){
        res_tot  = res
      }else{
        res_tot = rbind(res_tot, res)
      }
      a  = a+1
    }
  }
}

library(tidyverse)
library(ggplot2)
summary  = res_tot %>%
  group_by(prev,tau) %>%
  summarize(mean = mean(p), sd = sd(p), q97=quantile(p,0.975),q2=quantile(p,0.025))

summary$prev <- as.factor(summary$prev)
ggplot(summary)+
  geom_point(aes(x=tau, y=mean, colour = prev)) + 
  geom_ribbon(aes(x=tau, ymin=q2,ymax=q97, fill=prev), alpha=0.4 ) +
  theme_bw()