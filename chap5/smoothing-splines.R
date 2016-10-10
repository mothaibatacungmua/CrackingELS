library(MASS)
library(numDeriv)
library()
source('./utils.R')
source('./B-splines.R')


#####################################################
# Note: I haven't researched about the computation of smoothing spline yet.
# The below code just demonstrates in the ESL book but it's *very very* slow.
# The ESL book didn't mention to the caculation of the penalty term in the
# penalized residual sum of squares.

create_basis_func <- function(x, idx){
  basis = function(input){
    return(calc_b_spline(input, x, 4, 4, idx))
  }
  
  return(basis)
}

create_product_second_derv <- function(x, idx_i, idx_j){
  basis_i = create_basis_func(x, idx_i)
  basis_j = create_basis_func(x, idx_j)
  
  pen_func = function(input){
    ret = c()
    for(i in 1:length(input)){
      ret = c(ret, hessian(basis_i, input[i]) * hessian(basis_j, input[i]))
    }
    return(ret)
  }
  
  return(pen_func)
}

calc_pen_mat <- function(x){
  N = length(x)
  num_basis = N + 4
  pen_mat = matrix(rep(0, num_basis*num_basis), nrow = num_basis)
  for(i in 1:num_basis){
    for(j in 1:num_basis){
      pen_i_j = create_product_second_derv(x, i, j)
      pen_mat[i,j] = integrate(pen_i_j, min(x), max(x))[["value"]]
    }
  }
  
  return(pen_mat)
}

smoothing_spline_regression <- function(x, y, lambda){
  N = length(x)
  num_basis = N + 4
  
  b_mat = matrix(rep(0, N*num_basis), nrow = N)
  
  for(i in 1:num_basis){
    b_mat[,i] = b_spline(x, knots, order, order, i)
  }
  
  pen_mat = calc_pen_mat(x)
  reg_mat = ginv(crossprod(b_mat, b_mat) + lambda*pen_mat) %*% t(b_mat)
  smoother_mat = b_mat %*% reg_mat
  beta = reg_mat %*% y
  y_hat = smoother_mat %*% y
  
  return(list(beta=beta, y_hat=y_hat,smoother_mat=smoother_mat))
}
#####################################################

# At here, I use the smooth.spline function that was written in R to reproduce
# figure 5.6 in the book
data("bone")
View(bone)
male = which(bone$gender == "male")
female = which(bone$gender == "female")
plot(bone$age[male], bone$spnbmd[male], bg=1, col=1,pch=21, cex=0.5, xlab="Age", ylab="Relative Change in Spinal BMD")
fitted_male = smooth.spline(bone$age[male], bone$spnbmd[male], df=12)
lines(fitted_male$x, fitted_male$y, col=1, cex=2)

points(bone$age[female], bone$spnbmd[female], bg=2, col=2,pch=21, cex=0.5)
fitted_female = smooth.spline(bone$age[female], bone$spnbmd[female], df=12)
lines(fitted_female$x, fitted_female$y, col=2, cex=2)
