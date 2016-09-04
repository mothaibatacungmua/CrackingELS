# convert y
convert_class_vector <- function(y){
  f_y = as.factor(y)
  levels(f_y) = seq(1:nlevels(f_y))
  for(i in 1:nlevels(f_y)){
    y[which(f_y == i)] = i
  }
  return(y)
}

calc_prior_portions <- function(y){
  num_class = nlevels(as.factor(y))
  ret = c()
  for(i in 1:num_class){
    ret = c(ret, mean(y == i))
  }
  
  return(ret)
}

calc_vec_mu <- function(X, y){
  num_class = nlevels(as.factor(y))
  ret = matrix(rep(0, ncol(X)*num_class), nrow = num_class)
  
  for(i in 1:num_class){
    ret[i,] = colMeans(X[which(y == i),])
  }
  return(ret)
}

calc_common_s <- function(X, y, vec_mu){
  N = length(y)
  num_class = nlevels(as.factor(y))
  within = matrix(0, ncol(X), ncol(X))
  
  for(i in 1:num_class){
    X_t = X[which(y == i),]
    N_k = length(which(y == i))
    within = within + (N_k - 1)/(N - num_class)*var(X_t)
  }
  
  return(within)
}

# return a list of corvariance matrix for each class
calc_separate_s <- function(X, y, vec_mu){
  C = list()
  num_class = nlevels(as.factor(y))
  for(i in 1:num_class){
    X_t = X[which(y == i),]
    C[[i]] = var(X_t)
  }
  
  return(C)
}

# Page 113 in the ESL book, or see at: https://en.wikipedia.org/wiki/Whitening_transformation
whitenning_transform <- function(X, Sigma){
  ss = svd(Sigma)
  squared_d = diag(ss$d^(-1/2))
  trans = ss$u %*% squared_d
  return(list('X_w'= X %*% trans, 'mat_trans'=trans))
}

calc_between_class <- function(X, y){
  # counting N_k
  N = length(y)
  num_class = nlevels(as.factor(y))
  num_inst = c()
  for(i in 1:num_class){
    num_inst = c(num_inst, length(which(y==i)))
  }
  
  #mN = diag(num_inst/N, length(num_inst), length(num_inst))
  
  mat_mu = calc_vec_mu(X, y)
  common_mu = colMeans(X)
  
  between = matrix(0, ncol(X), ncol(X))
  
  for(i in 1:num_class){
    between = between + num_inst[i]/(N-1)* (mat_mu[i,] - common_mu) %*% t(mat_mu[i,] - common_mu)
  }
  return(between)
  
}

# do the same thing as calc_common_s function
calc_within_class <- function(X, y){
  N = length(y)
  num_class = nlevels(as.factor(y))
  #X_t = X
  #vec_mu = calc_vec_mu(X, y)
  within = matrix(0, ncol(X), ncol(X))
  
  for(i in 1:num_class){
    X_t = X[which(y == i),]
    N_k = length(which(y == i))
    within = within + (N_k - 1)/(N - num_class)*var(X_t)
  }
  
  return(within)
}