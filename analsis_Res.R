setwd("~/Dropbox/Group-testing/")
library(data.table)
source("compute_probas.R")

prev = 0.01
tau = 0.2
filename = "args[3]"testy_test
a = 0
N=15




for (alpha in seq(from=1, to=6, by=0.5) ){
  print(c(prev,tau))
  start_time <- Sys.time()
  res =rbindlist(lapply(1:300, function(it)
  {data.frame(p_unif = compute_p_tauprev_var(N, prev, tau, 0.1 * alpha, 0.05,
                                             B=10000, mode="uniform", mode_prev = "uniform"),
              p0 = compute_p_tauprev_var(N, prev, tau, 0.1* alpha, 0,
                                         B=10000, mode="uniform", mode_prev = "none"),
              p00 = compute_p_tauprev_var(N, prev, tau, 0.1* alpha, 0,
                                          B=10000, mode="none", mode_prev = "none"),
              p_mult = compute_p_tauprev_var(N, prev, tau, alpha, 0.05,
                                             B=10000, mode="multiplicative", mode_prev = "uniform"),
              p_mult = compute_p_tauprev_var(N, prev, tau, alpha, 0.05,
                                             B=10000, mode="multiplicative", mode_prev = "none"),
              prev = prev,
              tau=tau,
              var_mode = mode,
              var = alpha,
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
write.csv(res_tot, file = paste0("/scratch/users/cdonnat/Group_testing/experiments/", filename, ".csv" ))