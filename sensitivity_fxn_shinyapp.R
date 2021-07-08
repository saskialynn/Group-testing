# Compute sensitivity of the pooled PCR
# RV to simulate sensitivity 
library(fitdistrplus)
# Get all tests
tests <- data.table::fread("alltests_1mar24jun_v1.csv")

# Keep only valid positive ct values for first tests
tests %>% filter(result == "positive", firsttest==TRUE,!is.na(cttarget)) %>% pull(cttarget) -> cts

# Fit Weibull Distribution
fw <- fitdist(cts, "weibull")


## GENERATE Ct DATA
# Add variable of proportion Ct value >LoD to model, as surrogate for differences in viral load distribution in different populations

# Number of replicates to generate ecdf offset values for above.
set.seed(42)
x<-5000

# Create simulated vector of ct values with same shape as that from Weibull distribution fit to ct values from real first tests
# 5000 draws from weibull with shape and scale parameters drawn from real data
ct_fake_input <- rweibull(x,fw[[1]][1],fw[[1]][2]) 

# Create matrix of desired input parameters
# Change above.lod to % samples with ct value >LoD to reflect actual population of interest. Changing "lod" itself has no effect on model output.
lod = 35
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

### Include the probit  coefficients
# probit data input
probit_input <- read.csv("probit_zscores_cts_tissue_agnostic.csv")
probit_t<-subset(probit_input,probit_input$z_value<=1.96 & probit_input$z_value>=-1.96)
probit<-probit_t[,c(2,4,3)]
probit<-probit[order(probit$z_value,probit$ct_value),]
z_scores<-as.numeric(unlist(distinct(probit,z_value)))

# Number of replicates for model
set.seed(42)
# n <- 10000

probit.mode<-c("base","dsa.lower","dsa.upper","psa") 
probit.mode.index<-1 # 1 = no variation, 2, = LLN, 3 = ULN, 4 = probabilistic
probit.z.indices<-c(488,1,length(z_scores)) # 488 is a z score of 0 (base case) in the z index vector
dilution.vary.index<-1 # 1 = no variation, 2 = probabilistic


sens_fxn <- function(N, above.llod =35, B = 1000, ct_dat = ct_fake_input){
  above = which(mat2[,2] == above.llod/100 )
  set.seed(123)
  ct_set <- pmax(5,subset(ct_dat + mat2[above,3], ct_dat + mat2[above,3]<45)) 
  # (1) generate additional set of Ct values with above = 5-30% of values >LoD
  ## subset to only valid Ct values (<45)
  # (2) ct_set is the parallel maximum of the C_t values from (1) vs. 5 
  
  l <- sapply(1:N, function(positives){
    #threshold.ct <- c(sample(ct_set, 10000, replace=T)) # sample ct_set with replacement
    #dat <- matrix(sample(threshold.ct, positives * B, replace=T), nrow=positives) 
    dat <- matrix(sample(ct_set, positives * B, replace=T), nrow=positives) 
    # sample data uniformly at random 
    # positives*B samples of positives, rearrange into a matrix 
    each.conc = data.frame("x"= round(-log2(colSums(2^-dat)/N) + 
                                        ifelse(dilution.vary.index==1,0, 
                                               rnorm(mean=0,sd=1.1,n=ncol(dat))), 1))
    # sd of 1.1 reflects confidence interval for deviation from perfect log2 dilution in assays
    # calculation dilution based on number of positives (colSum) in total pool size
    # z.index= probit.z.indices[probit.mode.index]
    # probit score: probability of false positive or false negative test
    p = merge(probit %>% filter(z_value == 0), each.conc, by.x = "ct_value", by.y = "x", 
              all.y = TRUE, all.x = FALSE)
    # filter to z_value == 0 to get the median/mean sensitivity
    
    #return(mean(probit[1+(z.index-1)*571+each.conc*10-(lod-35.9)*10,2]))
    cbind(N, positives, mean(p$probit_probability_detection), 
          median(p$probit_probability_detection), 
          sd(p$probit_probability_detection), 
          as.numeric(quantile(p$probit_probability_detection, 0.025)), 
          as.numeric(quantile(p$probit_probability_detection, 0.975)), 
          mean(each.conc$x))
    # returns probability that a test is positive given k positives in sample
  })
  tl <- t(l) # transpose the matrix
  colnames(tl) <- c("pool.size", "positives", "mean", "median", 
                    "sd_mean", "q025", "q975", "ct_avg")
  return(tl)
} 
