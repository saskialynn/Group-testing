# Compute sensitivity of the pooled PCR
# RV to simulate sensitivity 
sens_fxn <- function(N, B = 1000, above = 5, ct_dat = ct_fake_input){ 
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