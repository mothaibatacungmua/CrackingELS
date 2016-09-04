gen_data <- function(N1, N2, N3){
  Sigma = matrix(c(1, 0, 0, 2), nrow=2)
  mu1 = c(7, 5) 
  mu2 = c(2, 2)
  mu3 = c(4, 7)
  
  set.seed(1)
  X_c1 = mvrnorm(N1, mu1, Sigma)
  X_c2 = mvrnorm(N2, mu2, Sigma)
  X_c3 = mvrnorm(N3, mu3, Sigma)
  X = rbind(X_c1, X_c2, X_c3)
  
  y_c1 = rep(0, N1)
  y_c2 = rep(1, N2)
  y_c3 = rep(2, N3)
  y = c(y_c1, y_c2, y_c3)
  
  return(list('X'=X, 'y'=y))
}

gen_train_data <- function(){
  return(gen_data(200, 200, 200))
}

gen_test_data <- function(){
  return(gen_data(50, 50, 50))
}