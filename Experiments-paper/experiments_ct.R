
library(tidyverse)
library(ggplot2)
library(fitdistrplus)
library(DescTools)
library(data.table)
library(dplyr)
library(reshape2)
library(RColorBrewer)


tests <- data.table::fread("~/Dropbox/Group-testing/alltests_1mar24jun_v1.csv")

# Keep only valid positive ct values for first tests
tests %>% filter(result == "positive", firsttest==TRUE,!is.na(cttarget)) %>% pull(cttarget) -> cts
mean(tests$cttarget, na.rm=TRUE)

plotdist(cts, histo = TRUE, demp = TRUE, breaks=18)
descdist(cts, boot = 1000)

ct_fake_input <- rweibull(x,fw[[1]][1],fw[[1]][2]) 

# Create matrix of desired input parameters
# Change above.lod to % samples with ct value >LoD to reflect actual population of interest. Changing "lod" itself has no effect on model output.
lod <- 35
above.lod <- seq(0.05,0.3,0.05)
translation_vector <- seq(-10,15,0.01)
mat<-matrix(ncol=3,nrow=length(translation_vector))

# Loop 1: fix LOD; for known shift, what percent of samples are above LoD (a)? 
for(v in 1:length(translation_vector)){
  i=1 
  fn<-ecdf(subset(pmax(5,ct_fake_input+translation_vector[v]),
                  (pmax(5,ct_fake_input+translation_vector[v]))<45)) ## selects adequate values within the sample
  a<-1-fn(lod) ## percent of cvalues in the distribution above the LoD 
  mat[(v-1)+i,1]<-lod
  mat[(v-1)+i,2]<-a
  mat[(v-1)+i,3]<-translation_vector[v]
}

# Shift ct values for each LoD and %above lod
mat2<-matrix(ncol=3,nrow=length(above.lod))

# Loop 2: fix LOD; for known a, how much should you shift? 
for(j in 1:length(above.lod)){
  tmp<-subset(mat,mat[,1]==lod)
  u<-tmp[which.min(abs(above.lod[j]-tmp[,2])),3]
  mat2[j,1]<-lod[i]
  mat2[j,2]<-above.lod[j]
  mat2[j,3]<-u
}


uu = c()
ct_dat <- rweibull(x,fw[[1]][1],fw[[1]][2])

ct_set <- pmax(5,subset(cts + mat2[above,3], cts  + mat2[above,3]<45)) 
res <- c()
for (positives in 1:50){
  res<- rbind(res,
             data.frame("type" = positives, 
                        "x"=sapply(1:50000, function(x){
    min(sample(cts, positives, replace=T))
  } )))
}

res2 = res %>% group_by(type) %>% summarise_all(mean)
print(res2)

ggplot(res %>% filter(type %in% c(1, 2,3,5, 10, 20))) +
  geom_density(aes(x=x, fill=as.factor(type)), alpha =0.5, adjust = 2) +theme_bw() +
  theme(legend.text=element_text(size=12))+labs(fill="# Positive \n Samples", size=12) +
  scale_fill_brewer(palette="Dark2")+
  ylab("Probability Density for the \n Minimum Ct in k positive samples")+
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=90, hjust=1))

  

uuu = res %>% group_by(type) %>% summarise_all(mean)

ggplot(res %>% filter(type %in% c("Individual Testing","N=1", "N=2", "N=3", "N=5", "N=10"))) +
  geom_histogram(aes(x=x, fill=type), alpha =0.5, binwidth = 2) +theme_bw()

uu <- c()
B=10000
for (N in c(2:20, 30, 40, 50, 60, 70, 80, 100)){
  ct_set <- pmax(5,subset(cts, cts<45)) 
  dat <- matrix(sample(ct_set, positives * B, replace=T), nrow=1 )
  # sample data uniformly at random 
  # n samples of positives, rearrange into a matrix of positive rows, n columns
  each.conc = data.frame("x"= sample(ct_set, B, replace = TRUE), type=paste0("Individual Testing"),
                         "Pool Size" = paste0("N=",N)); 
  for (positives in 2:min(50, N)){
    print(positives)
    # mean(sapply(ct_set, function(x){ uu$m[which(uu$ct_value == round(x))]}))
    # (1) generate additional set of C_t values with above = 5-30% of values above each LoD; subset to only valid Ct values (<45)
    # (2) ct_set is the parallel maximum of the C_t values from (1) vs. 5 (why 5??)
    
    dat <- matrix(sample(ct_set, positives * B, replace=T), nrow=positives) 
    # sample data uniformly at random 
    # n samples of positives, rearrange into a matrix of positive rows, n columns
    each.conc = rbind(each.conc, data.frame("x"=  -log2(colSums(2^-dat)/N), type=paste0("K=", positives),
                                            "Pool Size"= paste0("N=",N)))
    }
  
  each.conc$type = factor(each.conc$type, levels = c("Individual Testing", sapply(1:N, function(x){paste0("K=",x)})) )
  #ggplot(data=each.conc, aes(x=x, fill = type)) +
  #  geom_histogram( alpha=0.5)+ theme_bw()
  
  uu =rbind(uu,
            each.conc)
  
}

uu %>% group_by(type, Pool.Size) %>% summarise_all(mean)
uu$Pool.Size= factor(uu$Pool.Size, levels  = sapply(c(2:20, 30, 40, 50, 60, 70, 80, 100), function(x){paste0("N=",x)}))
ggplot(uu %>% filter(type %in% c("Individual Testing","K=1", "K=2", "K=3", "K=5", "K=10"), Pool.Size %in% c("N=2","N=8", "N=20", "N=50") )) +
  geom_histogram(aes(x=x, y=..count../sum(..count..), fill=type), alpha =0.5, binwidth = 2) +theme_bw() +
  facet_wrap(.~Pool.Size) + 
  scale_fill_brewer(palette="Dark2")+ xlab("") +  ylab("Density") +
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=90, hjust=1))
uuu = uu %>% group_by(type, Pool.Size) %>% summarise_all(mean)
uuu["adding"] = sapply(uuu$Pool.Size, function(x){log2(as.numeric(str_split(x, "=")[[1]][2]))})
uuu["min"] = sapply(uuu$type, function(x){ifelse(x=="Individual Testing", as.numeric((res2 %>% filter(type==1))["x"]),
                                                 as.numeric((res2 %>% filter(type==as.numeric(str_split(x, "=")[[1]][2])))["x"]) )})
uuu["Positives"] = sapply(uuu$type, function(x){ifelse(x=="Individual Testing", 1,
                                                       as.numeric(str_split(x, "=")[[1]][2])) })
ggplot(uuu) +
  geom_line(aes(x=type, y = min + adding, group = Pool.Size))
  
print(uu)
plot(x=1:N, (-1/(1:uu))/log(2))
points(x=1:(N-1),  diff(uu$mu[2:(N+1)]), col="red")
points(x=1:(N-1),  diff(uu$med[2:(N+1)]), col="blue")

plot(x=1:N, (-1/(1:N))/log(2))
points(x=1:(N-1),  diff(uu$med[2:(N+1)]), col="red")

proba <- sapply(1:N, function(positives){
  dat0 = 2^(-sample(ct_set,  B, replace=T))
  dat <- matrix(sample(ct_set, positives * B, replace=T), nrow=positives) 
  each.conc =  colSums(2^-dat)
  hist(dat0/each.conc)
  
})


mean(1/(1+2^(sapply(1:100000, function(b){ sample(cts,1)- sample(cts,1)})-1)))

sample(cts)
