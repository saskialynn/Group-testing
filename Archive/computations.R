library(tidyverse)
library(ggplot2)
library(fitdistrplus, lib.loc="~/R_libs")
library(DescTools, lib.loc="~/R_libs")
library(data.table, lib.loc="~/R_libs")
library(dplyr)
set.seed(42)

n<-5000
source("/scratch/users/cdonnat/Group_testing/Group-testing/compute_probas.R") 
source("/scratch/users/cdonnat/Group_testing/Group-testing/compute_proba_laws.R")
pool.max = 100
args = commandArgs(trailingOnly=TRUE)
prevalence = as.numeric(args[1])
tau = as.numeric(args[2])
pool = as.numeric(args[3])
mode = args[4]
mode_prev = args[5]
filename = args[6]
x<-50000

# probit data input
probit_input <- read.csv("probit_zscores_cts_tissue_agnostic.csv")
probit_t<-subset(probit_input,probit_input$z_value<=1.96 & probit_input$z_value>=-1.96)
probit<-probit_t[,c(2,4,3)]
probit<-probit[order(probit$z_value,probit$ct_value),]
z_scores<-as.numeric(unlist(distinct(probit,z_value)))

probit.mode<-c("base","dsa.lower","dsa.upper","psa") 
probit.mode.index<-1 # 1 = no variation, 2, = LLN, 3 = ULN, 4 = probabilistic
probit.z.indices<-c(488,1,length(z_scores)) # 488 is a z score of 0 (base case) in the z index vector
dilution.vary.index<-1 # 1 = no variation, 2 = probabilistic


tests <- data.table::fread("alltests_1mar24jun_v1.csv")
# Keep only valid positive ct values for first tests
tests %>% filter(result == "positive", firsttest==TRUE,!is.na(cttarget)) %>% pull(cttarget) -> cts


# Fit distributions
fw <- fitdist(cts, "weibull")
fno <- fitdist(cts, "norm")
fg <- fitdist(cts, "gamma")
fln <- fitdist(cts, "lnorm")

shape_w = summary(fw)[[1]][1]
scale_w = summary(fw)[[1]][2]


# Create simulated vector of ct values with same shape as that from Weibull distribution fit to ct values from real first tests
# 5000 draws from weibull with shape and scale parameters drawn from real data
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







experiment_loop <- function(above, ct_dat, pool.size, prev, tau){
  ct_set <- pmax(5,subset(ct_dat + mat2[above,3], ct_dat + mat2[above,3]<45)) 
  # (1) generate additional set of C_t values with above = 5-30% of values above each LoD; subset to only valid Ct values (<45)
  # (2) ct_set is the parallel maximum of the C_t values from (1) vs. 5 (why 5??)
  threshold.ct <- c(sample(ct_set, n, replace=T)) # sample ct_set with replacement # prevalence ("proportion positive tests")
  
  A = rbindlist(lapply(c(1, 1.1, 1.2, 1.5, 2, 2.5, 3,4,5,6,7,10), function(tau_relative_var) {
 print(c(pool.size, prev, tau, tau_relative_var, mode, mode_prev)) 

 
 sim_probs <- compute_probas_tauprev_var(pool.size, prev, tau, alpha=tau_relative_var, alpha_prev=0.05,
                                                     B=50000, mode =mode,
                                                     mode_prev=mode_prev)
 rbindlist(lapply(0:pool.size, function(positives) { # number of positives
#    print(positives)
          if (positives == 0) {
            # prevalence_corr <-  1-((1-prevalence)*(1-(prevalence*tau))^(p-1))
            data.frame(limit=lod, pool=pool.size, pos=positives, prevalence=prevalence, 
                       above.llod = above.lod[above], 
                       concentration=0, 
                       mode = mode,
                       mode_prev =mode_prev,
                       tau = tau,
                       prevalence = prev,
                       taus_tilde = sim_probs$tau[1],
                       tau_relative_var = tau_relative_var,
                       n_eff = sim_probs$n_eff[1],
                       N_eff = sim_probs$N_eff[1],
                       k_eff = sim_probs$k_eff[1],
                       probability = as.numeric(sim_probs[which(sim_probs$n == positives), "p"]),
                       #probability_subset1 =  sim_subset_probs[which(sim_subset_probs$n == 0), "p"],
                       probability_null = dbinom(positives, pool.size , prevalence), # uncorrected prevalence
                       probability_corr = dbinom(sim_probs$k_eff[1],ceiling(sim_probs$n_eff[1]), sim_probs$pi_eff[1]),
                       prevalence_corr = sim_probs$p[1],
                       pi = sim_probs$pi[1],
                       tau = sim_probs$tau[1],
                       probability_null = sim_probs$prob_null_eff[positives+1], # prevalence corrected for correlated individuals
                       random=0, 
                       z.index=0,
                       call.each.conc=FALSE, tests=1, tn=1,tp=0,fn=0,fp=0)
          } 
          else {
            
  
            dat <- matrix(sample(threshold.ct, positives * n, replace=T), nrow=positives) 
            # sample data uniformly at random 
            # n samples of positives, rearrange into a matrix of positive rows, n columns
            each.conc = -log2(colSums(2^-dat)/pool.size)+ifelse(dilution.vary.index==1,0,
                                                        rnorm(mean=0,sd=1.1,n=ncol(dat))) 
            # sd of 1.1 reflects confidence interval for deviation from perfect log2 dilution in assays
            # calculation dilution based on number of positives (colSum) in total pool size (p)
            z.index= probit.z.indices[probit.mode.index]
            
            tt = data.frame(
              limit=lod, pool=pool.size, pos=positives, prevalence=prevalence, 
              above.llod = above.lod[above], 
              concentration=0, 
              mode = mode,
              mode_prev =mode_prev,
              tau = tau,
              prevalence = prev,
              taus_tilde = sim_probs$tau[1],
              tau_relative_var = tau_relative_var,
              n_eff = sim_probs$n_eff[positives + 1],
              N_eff = sim_probs$N_eff[positives + 1],
              k_eff = sim_probs$k_eff[positives + 1],
              probability = as.numeric(sim_probs[which(sim_probs$n == positives), "p"]),
              #probability_subset1 =  sim_subset_probs[which(sim_subset_probs$n == 0), "p"],
              probability_null = dbinom(positives,pool.size , prevalence), # uncorrected prevalence
              probability_corr = dbinom(sim_probs$k_eff[positives + 1],ceiling(sim_probs$n_eff[1]), sim_probs$pi_eff[1]),
              prevalence_corr = sim_probs$p[1],
              pi = sim_probs$pi[1],
              tau = sim_probs$tau[1],
              probability_null = sim_probs$prob_null_eff[positives+1], # prevalence corrected for correlated individuals
              #random=sample(n,n,replace = TRUE)/n, #
              z.index=ifelse(probit.mode.index<4,probit.z.indices[probit.mode.index],
                             sample(1:length(z_scores),n,replace=T))) %>%
              mutate(
                call.each.conc=mean(probit[1+(z.index-1)*571+each.conc*10-(lod-35.9)*10,2]),
                tests=1 + pool * (call.each.conc), # number of tests done (number positive pools + pool size)
                tn=0,
                tp=1 * (call.each.conc),
                fn=1 * (!call.each.conc),
                fp=0)
            return(tt)
          }
        }))
  }))

}

allfirst.poolct <- experiment_loop(3, ct_fake_input, pool, prevalence, tau) # just for 15% Ct > LoD 

group_by(allfirst.poolct, pool, prevalence, above.llod, limit, tau,tau_relative_var, mode, mode_prev) %>% 
  summarize(pos1=weighted.mean(pos, w=probability),
            n_eff = mean(n_eff),
            prevalence_corr = mean(prevalence_corr),
            pos_corr=weighted.mean(pos, w=probability_corr),
            pos_null=weighted.mean(pos, w=probability_null),
            total.tests=weighted.mean(tests, w=probability), 
            total.tests_corr=weighted.mean(tests, w=probability_corr), 
            total.tests_null=weighted.mean(tests, w=probability_null), 
            tests.per.sample=weighted.mean(tests, w=probability)/mean(pool), # calculate tests per sample
            tests.per.sample_corr=weighted.mean(tests, w=probability_corr)/mean(pool),
            tests.per.sample_null=weighted.mean(tests, w=probability_null)/mean(pool),
            tn1=weighted.mean(tn, w=probability), 
            tp1=weighted.mean(tp, w=probability), 
            fn1=weighted.mean(fn, w=probability), 
            fp1=weighted.mean(fp, w=probability),
            ppa = weighted.mean(tp/(fn + tp), w=probability, na.rm = TRUE),
            ppa_corr = weighted.mean(tp/(fn + tp), w=probability_null, na.rm = TRUE),
            tn_corr=weighted.mean(tn, w=probability_corr), 
            tp_corr=weighted.mean(tp, w=probability_corr), 
            fn_corr=weighted.mean(fn, w=probability_corr), 
            fp_corr=weighted.mean(fp, w=probability_corr),
            ppa_null = weighted.mean(tp/(fn + tp), w=probability_null, na.rm = TRUE),
            tn_null=weighted.mean(tn, w=probability_null), 
            tp_null=weighted.mean(tp, w=probability_null), 
            fn_null=weighted.mean(fn, w=probability_null), 
            fp_null=weighted.mean(fp, w=probability_null)) -> apw


write.csv(apw, file = paste0("/scratch/users/cdonnat/Group_testing/experiments/", filename, ".csv" ))
