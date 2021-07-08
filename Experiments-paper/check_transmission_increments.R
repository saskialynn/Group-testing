library(tidyverse)
compute_pi_threshold <- function(pi, k, N){
  1/pi * ( 1-(1-k/(N*(1-pi)))^(1/(N-1)))
}


compute_increment <- function(pi, tau, N){
  N * ( 1-pi) * (1-(1-tau*pi)^(N-1))
}




plus<- c()
for (k in 1:10){
  print(k)
  for (pi in seq(from=0.001, to=0.2, by=0.001)){
    for (N  in 2:100){
      plus<- rbind(plus, c(pi, N, k,  compute_pi_threshold(pi, k, N)))
    }
  }
}
plus <- data.frame(plus)
colnames(plus)<- c("Prevalence", "Pool Size", "Increment", "Associated tau")
plus$`Associated tau`[which(plus$`Associated tau` >1)] = NA

dose.labs <- sapply(1:10, function(x){paste0("Mean Increment of ", x)})
names(dose.labs) <- 1:10


ggplot(plus %>% filter(Increment < 5))+
  geom_tile(aes(x=`Pool Size`, y =Prevalence, fill =  `Associated tau`)) +
  theme_bw() + 
  scale_fill_gradientn(trans = "log10", colours = terrain.colors(20)) +
  theme(legend.text=element_text(size=16))+labs(colour="Average Probability", fill="Minimum Value of \u03c4", size=14) + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1)) + 
  ylab("Prevalence π")  + facet_wrap(.~Increment, nrow=1, labeller = labeller(Increment = dose.labs)) 
