library(tidyverse)
library(ggplot2)

it=0
files <- list.files(path="~/Downloads/prevs", pattern="*.csv", full.names=TRUE, recursive=FALSE)
for (filename in files){
  ff = read_csv(filename)
  if (it == 0){
    res = ff
  }else{
    res = rbind(res, ff)
  }
  it = it+1
}


res$var =  as.factor(res$var)

summary  = res %>%
  dplyr::group_by(prev,tau, var) %>%
  dplyr::summarise(bias = mean(p), 
            sd = sd(p), 
            q97=quantile(p,0.975),
            q2=quantile(p,0.025))
summary  = res %>%
  group_by(prev,tau, var) %>%
  summarise(bias = mean(p -prev), 
            sd = sd(p -prev), 
            q97=quantile(p-prev,0.975),
            q2=quantile(p-prev,0.025))ana

summary$prev <- as.factor(summary$prev)
summary$var <- as.factor(summary$var)
#summary$prev <- as.factor(summary$prev)

ggplot(summary %>% filter( prev %in% c(0.001, 0.01, 0.02, 0.05, 0.1, 0.2)))+
  geom_point(aes(x=tau, y=bias, colour = prev)) + 
  geom_ribbon(aes(x=tau, ymin=q2,ymax=q97, fill=prev, colour=prev), alpha=0.4 ) +
  facet_wrap(vars(var)) + 
  theme_bw()

