# OBJECTIVES / METHOD
#' to create correlation within groups of size >2:
#' 
## Step 1: 

library(MASS) # contains mvrnorm fxn

create_correlated_Cts <- function(type="block", N=4, rho=0.8, shape = 4.55, scale = 29.96){
  # 1. create vector of length = size(group) by sampling from multivariate normal dist (mvrnorm())
  # # sample from multivariate normal distribution
  if (type == "random"){
    d = -1
    while(d<1e-3){
      Sigma <- matrix(runif(N^2, 0.3,0.9), 
                      nrow = N, ncol = N)
      Sigma = 0.5 * (Sigma + t(Sigma))
      diag(Sigma)= rep(1,N)
      d = det(Sigma)
    }
    
  }else{
    if (type == "block"){ 
      Sigma <-  (1-rho)  * diag(rep(1,N))  +rho * matrix(1, N, N)
      }
    if (type == "custom"){
      Sigma <- matrix(c(1, .1, .2, .3,
                        .1, 1, .01 , .4, 
                         .2, .01, 1, .5,
                          .3, .4, .5, 1), 
                         nrow = 4, ncol = 4)
      N=4}
  }
  # vector of means of the variables
  mu <- rep(0, N) 
  mvnorm_dat <- mvrnorm(1, mu, Sigma)
  
  ## Step 2: 
  #' 2. to each entry of vector apply inverse normal (pnorm) to get p-values
  p_val <- pnorm(mvnorm_dat)
  
  ## Step 3: 
  #' 3. Map back to distribution of interest (Weibull)
  weibull_dat <- qweibull(p_val, shape = shape, scale = scale)
  
  return(weibull_dat)
}























# Objective: generate data positively correlated with real data

####### A simple example #########
# https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables

######### Generating Correlated Data for Pooled Sampling #########
#' Notes
#' dont want all individuals correlated
#' create group assignments (column 'z')
#' within group correlation
#' procedure to generate groups: poisson distribution (centered around 1 or 2) determines group size
#' loop over different groups, within group generate correlated Ct values
#' Then imagine group assignment z hidden, shuffle data, repeat group testing model
#' what is the specificity/sensitivity under uniform testing?
#' now add extra information about group assignment; pool per group, repeat model 


# this code works for correlating pairs of people
#' to create correlation within groups of size >2:
#' create vector of length = size(group) by sampling from multivariate normal dist (mvrnorm())
#' map vector back to distribution of interest:
#' to each entry of vector apply inverse normal: qnorm()
#' 
#' GOAL: 
#' 
# Import Ct value info
#' probit_input <- read.csv("probit_zscores_cts_tissue_agnostic.csv")
#' 
#' # limit to just ct_values and first 10,000 rows of data 
#' # using more data produces "error: vector memory exhausted (limit reached)"
#' probit_input <- probit_input[c(1:10000),2] 
#' 
#' n <- length(probit_input)
#' rho <- 0.2                                  # desired correlation (= cos(angle))
#' theta <- acos(rho)                                         # corresponding angle
#' cor_dat <- rnorm(n, 1, 0.5)                                    # new random data
#' # should this be Weibull dis as well? which parameters to choose?
#' 
#' X <- cbind(probit_input, cor_dat) #matrix
#' Xctr <- scale(X, center = TRUE, scale = FALSE)     # centered columns (mean = 0)
#' rm(X)                                       # remove unneeded data to save space
#' 
#' Id <- diag(n)                                                # identity matrix
#' Q <- qr.Q(qr(Xctr[ , 1, drop = FALSE]))     # QR decomp. of Xctr (original data)
#' P <- tcrossprod(Q)             # = Q Q'    # projection onto space defined by x1
#' x2o <- (Id-P) %*% Xctr[ , 2]                    # x2ctr made orthogonal to x1ctr 
#' Xc2 <- cbind(Xctr[ , 1], x2o)                                   # bind to matrix
#' Y <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))            # scale columns to length 1
#' 
#' x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]      # final vector of correlated data
#' cor(probit_input, x)                                   # check correlation = rho
#' 
#' 
#' #' TODO: need to figure out how to do this for all Ct values 
#' #' (not just first 10,000 rows)
#' #' https://datasquad.at.sites.carleton.edu/data/storage-design/dealing-with-a-vector-memory-exhausted-error-in-r/
#' #' https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html
#' 
#' # working on whole data set
#' library(parallel)                          # base library for parallel computing
#' # Import Ct value info
#' probit_input <- read.csv("probit_zscores_cts_tissue_agnostic.csv")
#' 
#' # limit to just ct_values and first 10,000 rows of data 
#' # using more data produces "error: vector memory exhausted (limit reached)"
#' probit_input <- as.data.frame(probit_input[ ,2])
#' 
#' n <- nrow(probit_input)
#' rho <- 0.2                                  # desired correlation (= cos(angle))
#' theta <- acos(rho)                                         # corresponding angle
#' cor_dat <- as.data.frame(rnorm(n, 1, 0.5))                         # new random data
#' # should this be Weibull dis as well? which parameters to choose?
#' X <- cbind(probit_input, cor_dat) # matrix
#' Xctr <- scale(X, center = TRUE, scale = FALSE)     # centered columns (mean = 0)
#' rm(X)                                       # remove unneeded data to save space
#' 
#' fx <- function(Xctr){
#'   n <- nrow(Xctr)
#'   Id <- diag(n)  
#'   Q <- qr.Q(qr(Xctr[ , 1, drop = FALSE]))     # QR decomp. of Xctr (original data)
#'   P <- tcrossprod(Q)             # = Q Q'    # projection onto space defined by x1
#'   x2o <- (Id-P) %*% Xctr[ , 2]                    # x2ctr made orthogonal to x1ctr 
#'   Xc2 <- cbind(Xctr[ , 1], x2o)                                   # bind to matrix
#'   Y <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))            # scale columns to length 1
#'   return(Y)
#' }
#' 
#' system.time(
#'   results <- mclapply(Xctr, fx, mc.cores = detectCores()) # doesn't work :(
#' )
#' 
#' x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]      # final vector of correlated data
#' 
#' 
#' cor(probit_input, x)                                   # check correlation = rho
#' 

















