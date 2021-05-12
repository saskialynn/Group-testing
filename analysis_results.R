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

library(reshape2)

res2 = melt(res[, 2:ncol(res)], id.vars = c("prev", "tau", "var",  "var_mode", "it"))
summary  = res2 %>%
  dplyr::group_by(prev, tau, var, var_mode, variable) %>%
  dplyr::summarise(bias = mean(value), 
            sd = sd(value), 
            q97=quantile(value,0.975, na.rm = TRUE),
            q2=quantile(value,0.025, na.rm = TRUE))


test= res2 %>% filter(variable == "p00", tau ==0.5 , prev== 0.05)
test2= res2 %>% filter(variable == "p0", tau ==0.5 , prev== 0.05)


summary  = res %>%
  group_by(prev,tau, var) %>%
  summarise(bias = mean(p -prev), 
            sd = sd(p -prev), 
            q97=quantile(p-prev,0.975),
            q2=quantile(p-prev,0.025))

summary$prev <- as.factor(summary$prev)
summary$var <- as.factor(summary$var)
#summary$prev <- as.factor(summary$prev)

ggplot(summary %>% filter(variable == "p00"))+
  geom_point(aes(x=tau, y=sd, colour = prev)) + 
  #geom_ribbon(aes(x=tau, ymin=q2,ymax=q97, fill=prev, colour=prev), alpha=0.4 ) +
  facet_wrap(vars(var)) + 
  theme_bw()

