radial_basis <- function(x, center){
  dist = sqrt(sum((x-center)^2))
  if(dist == 0){
    return(0)
  }
  
  return(dist^2 * log(dist))
}

gen_kernel_matrix <-function(X){
  nr= nrow(X)
  ker_mat = matrix(rep(0, nr*nr), nrow=nr)
  
  for(i in 1:nr){
    for(j in 1:nr){
      ker_mat[i,j] = radial_basis(X[i,], X[j,])
    }
  }
  
  return(ker_mat)
}

# Note: We can instead of lambda by degree of freedom of the smoother matrix S by coding more :D
thin_plate_regression <- function(X, y, lambda){
  #browser()
  ker_mat = gen_kernel_matrix(X)
  input_mat = cbind(X, ker_mat)
  input_mat = cbind(rep(1, nrow(X)), input_mat)
  
  ker_mat = cbind(rep(0, nrow(input_mat)),rep(0, nrow(input_mat)),rep(0, nrow(input_mat)), ker_mat)
  ker_mat = rbind(ker_mat, rep(0, ncol(input_mat)),rep(0, ncol(input_mat)), rep(0, ncol(input_mat)))
  
  beta = ginv(crossprod(input_mat, input_mat) + lambda*ker_mat) %*% t(input_mat) %*% y
  y_hat = input_mat %*% beta
  
  return(list(beta=beta, y_hat=y_hat))
}

library(ElemStatLearn)
library(rgl)
data("SAheart")
X = cbind(SAheart$age, SAheart$obesity)
y = SAheart$sbp

ret = thin_plate_regression(X, y, 0.0001)
plot3d(X[,1], X[,2],y, pch=16, highlight.3d=TRUE, main="3D Scatterplot")
plot3d(X[,1], X[,2],ret$y_hat, pch=16, highlight.3d=TRUE, main="3D Scatterplot")
