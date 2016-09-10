# At here, I write the Newton method for the logistic regression 
# for n-class (Exercise 4.4 in the ESL book)
library(MASS)
library(base)
source('utils.R')
source('gen-data.R')


calc_softmax_prob <- function(X, theta){
  softmax_prob = X %*% theta
  softmax_prob = exp(softmax_prob)
  
  #remove infinite
  softmax_prob = apply(softmax_prob, 1, function(x){
    x[which(is.infinite(x))] = .Machine$double.xmax
    return(x)
  })
  
  softmax_prob = t(softmax_prob)
  
  denom_vec = rowSums(softmax_prob)
  
  denom_vec[which(is.infinite(denom_vec))] = .Machine$double.xmax
  
 
  for(i in 1:length(denom_vec)){
    softmax_prob[i,] = softmax_prob[i,]/denom_vec[i]
  }
  
  
  return(softmax_prob)
}

convert_yy <- function(y){
  n = nlevels(as.factor(y))
  conv = matrix(rep(0, length(y)*n), nrow = length(y))
  
  for(i in 1:length(y)){
    conv[i, y[i]] = 1
  }
  
  return(conv)
}

cost_fn <- function(X, y, theta){
  softmax_prob = calc_softmax_prob(X, theta)
  ind = convert_yy(y)
  
  cost = 0
  for(i in 1:length(y)){
    # This stupid code prevent log(0) * 0 = NaN
    for(j in 1:length(ind[i,])){
      if(softmax_prob[i,j] != 0){
        cost = cost + ind[i,j] * log(softmax_prob[i,j])
      }
    }
  }
  
  return(-cost)
}

calc_grad <- function(X, y, theta){
  log_prob = calc_softmax_prob(X, theta)
  
  g = theta
  
  # Note: We don't need caculate the gradient of beta_K
  # that is the last column of theta
  for(i in 1:(ncol(theta) - 1)){
    p = log_prob[, i]
    p[which(y == i)] = 1 - p[which(y == i)]
    p[which(y != i)] = -p[which(y != i)]
    g[, i] = crossprod(X, p)
  }
  
  return(g)
}

calc_hessian <- function(X, theta){
  log_prob = calc_softmax_prob(X, theta)
  w = log_prob * (1 - log_prob)
  
  hess = list()
  for(i in 1:(ncol(theta)-1)){
    mat = matrix(rep(0,ncol(X)*ncol(X)), ncol=ncol(X))
    for(j in 1:nrow(X)){
      mat = mat + X[j,] %*% t(X[j,]) * w[j,i]
    }
    hess[[i]] = -ginv(mat)
  }
  
  hess[[ncol(theta)]] = diag(ncol(X))
  
  return(hess)
}

# Note: Newton method is very slow when X_train is large. 
# Instead, we can use the conjugate gradient method or LBFS to solve fast 
# the logistic regession
set.seed(1)
logis_reg <- function(X, y, max_iter=20){
  nl = nlevels(as.factor(y))
  
  theta = matrix(runif(ncol(X) * nl), nrow = ncol(X))
  theta[,nl] = rep(0, ncol(X))
  
  for(i in 1:max_iter){
    grad = calc_grad(X, y, theta)
    hess = calc_hessian(X, theta)
    # browser()
    # monitor cost:
    print(sprintf('Iteration:%d Cost:%0.5f', i, cost_fn(X, y, theta)))
    
    
    for(j in 1:nl){
      theta[,j] = theta[,j] - hess[[j]] %*% grad[,j]
    }
  }
  
  return(theta)
}

logis_pred <- function(X, theta){
  log_prob = calc_softmax_prob(X, theta)
  pred = rep(1, nrow(X))
  for(i in 1:nrow(X)){
    pred[i] = which.max(log_prob[i,])
  }
  return(pred)
}

wrapper <- function(X_train, y_train, X_test, y_test){
  theta = logis_reg(X_train, y_train)
  train_pred = logis_pred(X_train, theta)
  test_pred = logis_pred(X_test, theta)
  
  return(c(mean(train_pred == y_train), mean(test_pred == y_test)))
}


train_data = gen_train_data()
test_data = gen_test_data()
X_train = train_data$X
y_train = convert_class_vector(train_data$y)

X_test = test_data$X
y_test = convert_class_vector(test_data$y)

X_train = cbind(rep(1, nrow(X_train)), X_train)
X_test = cbind(rep(1, nrow(X_test)), X_test)


wrapper(X_train, y_train, X_test, y_test)

