prev=0.01
res = lapply(seq(from=0.001, to=0.6, by =0.001), function(tau){
  t(sapply(seq(from=2,to=50, by =1), function(N){
    c(tau, N, compute_corr(N, prev, tau))
  }))
  
})
res = do.call(rbind, res)
res = data.frame(res)
colnames(res) = c("tau", "N", "corr")




library(ggplot2)
ggplot(res, aes(x=N, y=tau, fill=corr))+
  scale_fill_gradient(name = expression(rho), low = "blue", high="red")+
  geom_tile() +
  ylab(expression(tau))+
  xlab("Pool size N")+
  scale_x_discrete(breaks = seq(from=5, to=50, by=5))+
  theme_bw() #+
#theme(axis.text.x=element_blank(),
#axis.ticks.x=element_blank())

library(tidyverse)
res$N <- as.factor(res$N)
ggplot(res %>% filter(N %in% c(5,10,15,20,30)), aes(x = tau, y=corr,colour=N ))+
  geom_line(size=1.2)+
  theme_bw() + ylab(expression(rho)) + xlab(expression(tau))


prev=0.01
res = lapply(seq(from=0.001, to=0.6, by =0.001), function(tau){
  t(sapply(seq(from=2,to=50, by =1), function(N){
    c(tau, N, compute_p(N, prev, tau))
  }))
  
})
res = do.call(rbind, res)
res = data.frame(res)
colnames(res) = c("tau", "N", "p")
res$N <- as.factor(res$N)
ggplot(res %>% filter(N %in% c(5,10,15,20,30,40,50)), aes(x = tau, y=p,colour=N ))+
  geom_line(size=1.2)+
  theme_bw() + xlab(expression(tau))

#### Kolmogorov distance:
res2 = lapply(seq(from=0.001, to=0.6, by =0.001), function(tau){
  t(sapply(seq(from=2,to=50, by =1), function(N){
    a = lapply(seq(from=0.0001, to=0.10, by=0.0001), function(prev){
      p1 = compute_probas(N, prev, tau)
      d = lapply(seq(from=1, to= N, by=1), function(n0){
        p2 = sapply(0:N, function(x){dbinom(x,n0, prev)})
        return(c(prev, tau, N, max(abs(p1$p - p2))))})
      d = do.call(rbind, d)
      return(d)
    })
    return(a)
  }))
})


N=10
res2 = lapply(seq(from=0.001, to=0.6, by =0.001), function(tau){
  do.call(rbind, lapply(seq(from=0.0001, to=0.10, by=0.001), function(prev){
    d = do.call(rbind,lapply(seq(from=2, to= 50, by=1), function(n0){
      p1 = compute_probas(n0, prev, tau)
      p2 = sapply(0:n0, function(x){dbinom(x,n0, prev)})
      return(c(prev, tau, n0, max(abs(p1$p - p2))))}))
  }))
})
res2 = do.call(rbind, res2)
res2 = data.frame(res2)
colnames(res2) = c("prev", "tau", "N", "distance")
ggplot(res2 %>% filter(N %in% c(2, 5,10,20,30)), aes(x=prev, y=tau, fill=distance))+
  scale_fill_gradient(name = "Distance", low = "blue", high="red")+
  geom_tile() +
  ylab(expression(tau))+
  xlab(expression(pi))+
  facet_wrap(~N, nrow=1)+
  #scale_x_discrete(breaks = seq(from=0.0001, to=50, by=0.0001))+
  theme_bw() 


res2$tau <- as.factor(res2$tau)
ggplot(res %>% filter(tau %in% c(0.001,0.01, 0.05, 0.1, 0.2,0.3,0.5,0.6)), 
       aes(x = prev, y=distance,colour=tau ))+
  geom_line(size=1.2)+
  theme_bw() + ylab(expression(rho)) + xlab(expression(tau))




res2 = lapply(seq(from=0.001, to=0.6, by =0.001), function(N){
  do.call(rbind, lapply(seq(from=0.0001, to=0.10, by=0.001), function(prev){
    d = do.call(rbind,lapply(seq(from=2, to= 50, by=1), function(n0){
      p1 = compute_probas(n0, prev, tau)
      p2 = sapply(0:n0, function(x){dbinom(x,n0, prev)})
      return(c(prev, tau, n0, max(abs(p1$p - p2))))}))
  }))
})


res3 =  lapply(seq(from=0.0001, to=0.10, by=0.001), function(prev){
  d = do.call(rbind,lapply(seq(from=2, to= 50, by=1), function(n0){
    p1 = compute_probas(n0, prev, tau)
    k = do.call(rbind, lapply(1:n0, function(ne){
      t(sapply(seq(from=0.0001, to=1.0, by=0.001), function(prev_eff){
        p2 = sapply(1:ne, function(x){dbinom(x,ne, prev_eff)})
        p2 = c(p2, rep(0, n0-ne))
        c(prev_eff, ne, max(abs(p1$p - p2)))
      }))}))
    ind = which.min(k[,3])
    return(c(prev, tau, n0, k[ind,]))
  }))
})
res3 = do.call(rbind,res3)
res3 = data.frame(res3)
colnames(res3) < c("prev", "tau", "N", "prev_e", "N_e", "d" )
res4 =  lapply(seq(from=0.0001, to=0.10, by=0.001), function(prev){
  t(sapply(seq(from=2, to= 50, by=1), function(n0){
    a = compute_effective_parameters2(n0, prev, tau)
    return(c(prev, tau, n0, a[1], a[2], a[3]))
  })
  )
})
res4 = do.call(rbind,res4)
res4 = data.frame(res4)
colnames(res4) <- c("prev", "tau", "N", "prev_e", "N_e", "d" )


res3$N =as.factor(res3$N )
ggplot(res3 %>% filter(N %in% c(5,10,15,20,50)), 
  aes(x = prev, y=N_e , colour=N))+
  geom_line(size=1.2)+
  theme_bw() 
ggplot(res3 %>% filter(N %in% c(5,10,15,20,50)), 
       aes(x = prev, y=prev_e , colour=N))+
  geom_line(size=1.2)+
  theme_bw() 

res4$N =as.factor(res4$N )
ggplot(res4 %>% filter(N %in% c(5,10,15,20,50)), 
       aes(x = prev, y=N_e , colour=N))+
  geom_line(size=1.2)+
  theme_bw() 

res4$N =as.factor(res4$N )
ggplot(res4 %>% filter(N %in% c(5,10,15,20,50)), 
       aes(x = prev, y=prev_e , colour=N))+
  geom_line(size=1.2)+
  theme_bw() 


ggplot(res %>% filter(tau %in% c(0.001,0.01, 0.05, 0.1, 0.2,0.3,0.5,0.6)), 
       aes(x = prev, y=distance,colour=tau ))+
  geom_line(size=1.2)+
  theme_bw() + ylab(expression(rho)) + xlab(expression(tau))


