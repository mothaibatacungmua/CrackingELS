library(MASS)
library(base)

source('perceptron.R')

dat = gen_two_class_data()
X = dat$X
y = dat$y

plot(X[,2], X[,3], col=as.numeric(y)+2, xlab = "X1", ylab = "X2", bg=as.numeric(y)+2, pch=21)

# This below algorithm was described at page 131 in the ESL book

cost_fn <- function(X, y, beta){
  y_hat = X %*% beta * y #elementwise matrix mulplication
  miss = which(y_hat < 0)
  
  cost = -sum(y_hat)
  
  return(cost)
}

grad <- function(X ,y, beta){
  y_hat = X %*% beta * y
  miss = which(y_hat < 0)
  
  # gradient for beta1
  grad_beta1 = -sum(X[miss,]*y[miss])
  
  # gradient for beta_0
  grad_beta0 = -sum(y[miss])
  
  return(c(grad_beta1, grad_beta0))
}


# this function is written according to the exercise 4.6 in the ESL, you can use cost_fn function 
# and grad function together with a learning rate to accomplish the original version of the perceptron algorithm
# Note: You should guaranteed that the data is separatable if not the algorithm won't be converge
perceptron <- function(X, y){
  #random init beta
  beta = rnorm(ncol(X),sd=6)
  y_hat = X %*% beta * y 
  miss = which(y_hat < 0)
  
  while(length(miss) > 0){
    beta = beta + y[miss[1]]*X[miss[1],]
    y_hat = X %*% beta * y 
    miss = which(y_hat < 0)
  }
  
  return(beta)
}

# you can see that the problem will have multiple solutions with the perceptron algorithms
set.seed(1)
beta = perceptron(X, y)
abline(-beta[1]/beta[3], -beta[2]/beta[3],col="black")

set.seed(10)
beta = perceptron(X, y)
abline(-beta[1]/beta[3], -beta[2]/beta[3],col="red")

set.seed(50)
beta = perceptron(X, y)
abline(-beta[1]/beta[3], -beta[2]/beta[3],col="blue")
