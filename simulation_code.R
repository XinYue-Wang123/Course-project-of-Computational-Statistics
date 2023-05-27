
############ Simulation Example 1: (In)Consistency with 3 Variables ############
rm(list=ls())
library(Matrix)
library(MASS)
library(glmnet)

####=== S1.1 data generating process ===####

mydata1 <- function(n, miu, sigma, beta){
  x.1    <- rnorm(n, miu, sigma)
  x.2    <- rnorm(n, miu, sigma)
  x.3    <- 2/3*x.1 + 2/3*x.2 + 1/3*rnorm(n,miu,sigma)
  x      <- cbind(x.1, x.2, x.3)
  error  <- rnorm(n, miu, sigma)
  y      <- x.1*beta[1] + x.2*beta[2]+ error
  data   <- data.frame(y, x)
  return(data)
}


####=== S1.2 setting a: ¦Â=(2,3) ===####
set.seed(123)
traindata.a <- mydata1(n = 1000, miu = 0, sigma = 1, beta = c(2, 3))
x           <- cbind(traindata.a$x.1, traindata.a$x.2, traindata.a$x.3)

#=== verify Strong Irrepresentable Condition (fails for setting a) ===#
X1 <- cbind(traindata.a$x.1, traindata.a$x.2)
X2 <- traindata.a$x.3
solve(t(X1) %*% X1) %*% t(X1) %*% X2

#=== the Lasso paths for setting a ===#
lasso.sim.a <- glmnet(x, traindata.a$y, family = 'gaussian',
                      alpha = 1, nlambda = 100,
                      standardize = T, intercept = F)
lasso.sim.a
plot(lasso.sim.a, lwd = 2)

cv.sim.a <- cv.glmnet(x, traindata.a$y, family = 'gaussian',
                      nfolds = 10, alpha = 1, nlambda = 100,
                      standardize = T, intercept = F)
cv.sim.a
plot(cv.sim.a, lwd = 2)


####=== S1.3 setting b: ¦Â=(-2,3) ===####
set.seed(123)
traindata.b <- mydata1(n = 1000, miu = 0, sigma = 1, beta = c(-2, 3))
x           <- cbind(traindata.b$x.1, traindata.b$x.2, traindata.b$x.3)

#=== verify Strong Irrepresentable Condition (holds for setting b) ===#
X1 <- cbind(traindata.b$x.1, traindata.b$x.2)
X2 <- traindata.b$x.3
solve(t(X1) %*% X1) %*% t(X1) %*% X2

#=== the Lasso paths for setting 2 ===#
lasso.sim.b <- glmnet(x, traindata.b$y, family = 'gaussian',
                      alpha = 1, nlambda = 100,
                      standardize = T, intercept = F)
lasso.sim.b
plot(lasso.sim.b, lwd = 2)

cv.sim.b <- cv.glmnet(x, traindata.b$y, family = 'gaussian',
                      nfolds = 10, alpha = 1, nlambda = 100,
                      standardize = T, intercept = F)
cv.sim.b
plot(cv.sim.b, lwd = 2)



############ Simulation Example 2: Quantitative Evaluation of Impact ############
library(Matrix)
library(MASS)
library(mvtnorm)        # package for X distribution
library(glmnet)         # package for the Lasso
library(matrixsampling) # package for Wishart distribution


###=== data generating process ===###

mydata2 <- function(n, p, q, beta, S){
  miu   <- 0
  sigma <- sqrt(0.1)
  error <- rnorm(n, miu, sigma) 
  X     <- rmvnorm(n, rep(miu,p), S)
  Y     <- X %*% beta + error
  data  <- data.frame(Y, X)
  return(data)
}

n    <- 100    # number of observations
p    <- 32     # number of variables
q    <- 5      # number of non-zero beta
beta <- c(7, 4, 2, 1, 1, rep(0,p-q))


###=== Impact of Strong Irrepresentable Condition on Model Selection ===###

correct    <- matrix(NA, 1000, 1)
percentage <- matrix(NA, 100, 1) # creat store for y-axis
eta        <- matrix(NA, 100, 1) # creat store for x-axis

for (i in 1:100){ # generate 100 design of X #
  set.seed(123)
  S <- rwishart(100, nu = 13, Sigma = diag(p)) [, , i]

  # verify Strong Irrepresentable Condition by ¦Ç #
  sim.data <- mydata2(n, p, q, beta, S)
  
  X1  <- as.matrix( sim.data[ , 2:(q+1)] )     # X with non-zero ¦Â
  X2  <- as.matrix( sim.data[ , (q+2):(p+1)] ) # X with zero ¦Â
  C11 <- 1/n * t(X1) %*% X1
  C21 <- 1/n * t(X2) %*% X1
  beta.sign <- as.matrix(sign(beta[1:q]))
  
  eta[i,1] <- 1 - norm(C21 %*% solve(C11) %*% beta.sign, type = 'I')
  
  # compute percentage of generating matched models #
  set.seed(123)
  for (T in 1:1000){ # generate 1000 simulations for each design #
    data <- mydata2(n, p, q, beta, S)
    
    # calculate the Lasso path #
    lasso.sim <- glmnet(as.matrix(data[,2:(p+1)]), data$Y,
                        alpha = 1, intercept = F)
    cv.sim <- cv.glmnet(as.matrix(data[,2:(p+1)]), data$Y,
                        alpha = 1, intercept = F)
    beta.sim <- predict(lasso.sim, type = 'coefficients',
                        s = cv.sim$lambda.min) [2:(q+1), ]
    
    # examine sign of estimate #
    ifelse(sum(sign(beta.sim)) == q, correct[T,1] <- 1, correct[T,1] <- 0)
    }
  percentage[i,1] <- sum(correct)/1000
}

# plot #
plot(eta, type = 'p', pch = 20, xlab=' ', ylab='¦Ç¡Þ')
abline(h = 0, lty = 2)
plot(eta,percentage, type = 'p', pch = 20, 
     xlab='¦Ç¡Þ', ylab='Percentage of Lasso Selecting the Correct Model')
abline(v = 0, lty = 2)

