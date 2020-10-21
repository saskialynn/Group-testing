
# Objective: generate data positively correlated with real data
# https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables

####### A simple example #########
n     <- 20                    # length of vector
rho   <- 0.6                   # desired correlation = cos(angle)
theta <- acos(rho)             # corresponding angle
x1    <- rnorm(n, 1, 1)        # fixed given data
x2    <- rnorm(n, 2, 0.5)      # new random data
X     <- cbind(x1, x2)         # matrix
Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)

Id   <- diag(n)                               # identity matrix
Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q (original data)
P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr (generated correlated data)
Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
cor(x1, x)                                    # check correlation = rho


######### Generating Correlated Data for Pooled Sampling #########
# Import Ct value info
probit_input <- read.csv("probit_zscores_cts_tissue_agnostic.csv")
probit_input <- probit_input[c(1:10000),2] # limit to just ct_values
n <- length(probit_input)
rho <- 0.2 # desired correlation (= cos(angle))
theta <- acos(rho) # corresponding angle
cor_dat <- rnorm(n, 1, 0.5) # new random data
# should this be Weibull dis as well? which parameters to choose?

X <- cbind(probit_input, cor_dat) #matrix
Xctr <- scale(X, center = TRUE, scale = FALSE) # centered columns (mean = 0)
rm(X)
rm(probit_input) # remove unneeded data to save space

Id <- diag(n) # identity matrix
Q <- qr.Q(qr(Xctr[ , 1, drop = FALSE])) 
# compute QR decomposition of Xctr (select column of original data)
P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr (generated correlated data)
Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
cor(x1, x)                                    # check correlation = rho



