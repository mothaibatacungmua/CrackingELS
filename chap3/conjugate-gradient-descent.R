# I just write this function for reference. 
gram_schmidt_conjugation <- function(A){
  n = ncol(A)
  U = diag(n)
  
  conju = matrix(rep(0, n*n), ncol = n)
  browser()
  conju[,1] = U[,1]
  for(i in 2:n){
    ex = rep(1, i-1)
    c_m = as.matrix(conju[,1:(i-1)])
    beta = -(t(c_m) %*% A %*% U[,i]) / 
           (t(c_m)  %*% A %*% c_m %*% ex)
    conju[,i] = U[,i] + c_m %*% beta
  }
  
  return(conju)
}

# Note: A must be a positive definite matrix
conjugate_gradient <- function(A, b){
  n = ncol(A)
  xi = runif(n)
  di = r_prev = b - A %*% xi
  #browser()
  for(i in 2:n){
    alpha = as.numeric(crossprod(r_prev, r_prev) / crossprod(di, A %*% di))
    xi = xi + alpha * di
    r_next = r_prev - alpha * A %*% di
    beta = as.numeric(crossprod(r_next, r_next) / crossprod(r_prev, r_prev))
    di = r_next + beta * di
    r_prev = r_next
  }
  
  return(xi)
}

# Note:
# Setting b = t(X) %*% y, A = t(X) %*% X, x_0 = 0. We can derive 
# PLS from the conjugate gradient sequence:
# x %*% d_i = z_m and alpha_i = <z_m, y> / <z_m, z_m>
# (Exercies 3.18 in the ESL book)

