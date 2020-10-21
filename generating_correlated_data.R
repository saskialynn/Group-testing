
# Objective: generate data positively correlated with real data

####### A simple example #########
# https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables

######### Generating Correlated Data for Pooled Sampling #########
# Import Ct value info
probit_input <- read.csv("probit_zscores_cts_tissue_agnostic.csv")

# limit to just ct_values and first 10,000 rows of data 
# using more data produces "error: vector memory exhausted (limit reached)"
probit_input <- probit_input[c(1:10000),2] 

n <- length(probit_input)
rho <- 0.2                                  # desired correlation (= cos(angle))
theta <- acos(rho)                                         # corresponding angle
cor_dat <- rnorm(n, 1, 0.5)                                    # new random data
# should this be Weibull dis as well? which parameters to choose?

X <- cbind(probit_input, cor_dat) #matrix
Xctr <- scale(X, center = TRUE, scale = FALSE)     # centered columns (mean = 0)
rm(X)                                       # remove unneeded data to save space

Id <- diag(n)                                                # identity matrix
Q <- qr.Q(qr(Xctr[ , 1, drop = FALSE]))     # QR decomp. of Xctr (original data)
P <- tcrossprod(Q)             # = Q Q'    # projection onto space defined by x1
x2o <- (Id-P) %*% Xctr[ , 2]                    # x2ctr made orthogonal to x1ctr 
Xc2 <- cbind(Xctr[ , 1], x2o)                                   # bind to matrix
Y <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))            # scale columns to length 1

x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]      # final vector of correlated data
cor(probit_input, x)                                   # check correlation = rho


#' TODO: need to figure out how to do this for all Ct values 
#' (not just first 10,000 rows)
#' https://datasquad.at.sites.carleton.edu/data/storage-design/dealing-with-a-vector-memory-exhausted-error-in-r/
#' https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html

# working on whole data set
library(parallel)                          # base library for parallel computing
# Import Ct value info
probit_input <- read.csv("probit_zscores_cts_tissue_agnostic.csv")

# limit to just ct_values and first 10,000 rows of data 
# using more data produces "error: vector memory exhausted (limit reached)"
probit_input <- as.data.frame(probit_input[ ,2])

n <- nrow(probit_input)
rho <- 0.2                                  # desired correlation (= cos(angle))
theta <- acos(rho)                                         # corresponding angle
cor_dat <- as.data.frame(rnorm(n, 1, 0.5))                         # new random data
# should this be Weibull dis as well? which parameters to choose?
X <- cbind(probit_input, cor_dat) # matrix
Xctr <- scale(X, center = TRUE, scale = FALSE)     # centered columns (mean = 0)
rm(X)                                       # remove unneeded data to save space

fx <- function(Xctr){
  n <- nrow(Xctr)
  Id <- diag(n)  
  Q <- qr.Q(qr(Xctr[ , 1, drop = FALSE]))     # QR decomp. of Xctr (original data)
  P <- tcrossprod(Q)             # = Q Q'    # projection onto space defined by x1
  x2o <- (Id-P) %*% Xctr[ , 2]                    # x2ctr made orthogonal to x1ctr 
  Xc2 <- cbind(Xctr[ , 1], x2o)                                   # bind to matrix
  Y <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))            # scale columns to length 1
  return(Y)
}

system.time(
  results <- mclapply(Xctr, fx, mc.cores = detectCores()) # doesn't work :(
)

x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]      # final vector of correlated data


cor(probit_input, x)                                   # check correlation = rho


















