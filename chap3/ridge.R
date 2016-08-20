library(MASS)

ridge <- function(X, y, lambda){
  A = crossprod(X, X)
  n = nrow(A)
  return((ginv(A + lambda*diag(n)) %*% crossprod(X,y)))  
}

# I want to solve the ridge by the BFGS method, you can dismiss it