library(MASS)

gen_data <- function(x){
  n = length(x)
  return(4*sin(-6*x+2) + rnorm(n, 0, 2))
}

x = seq(from=0, to=1, length.out = 50)
y = gen_data(x)
real_y = 4*sin(-6*x+2)
plot(x, real_y, xlab = "X", ylab = "Y", type="l", ylim = c(-6,6), cex=1.5)
points(x,y)


calc_cubic_spline_basis <- function(x, order, knot){
  if(order == 1){
    return(rep(1, length(x)))
  }
  
  #browser()
  
  if(is.nan(knot)){
    return(x^(order-1))
  }
  
  z = rep(0, length(x))
  v = (x - knot)^(order - 1)
  ind = which(v > z)
  z[ind] = v[ind]
  return(z)
  
}

calc_cubic_spline_mat <- function(x, order, knots){
  n = length(x)
  imat = matrix(rep(1,n), ncol = 1)
  
  #browser()
  for(i in 2:order){
    imat = cbind(imat, calc_cubic_spline_basis(x, i, NaN))
  }
  
  for(i in 1:length(knots)){
    imat = cbind(imat, calc_cubic_spline_basis(x, order, knots[i]))
  }
  
  return(imat)
}

cubic_spline_regression <- function(x, y, order, nknots){
  knots = seq(min(x), max(x), length.out = (nknots+2))
  knots = knots[2:(length(knots)-1)]
  
  imat = calc_cubic_spline_mat(x, order, knots)
  beta = ginv(crossprod(imat, imat)) %*% t(imat) %*% y
  y_hat = imat %*% beta
  
  return(list('knots'=knots, 'beta'=beta, 'y_hat'=y_hat))
}

ret = cubic_spline_regression(x, y, order=4, nknots = 2)
lines(x, ret$y_hat, col="red")
