library(MASS)

source('./utils.R')
dat = gen_data()
x = dat$x
y = dat$y
real_y = dat$real_y

plot(x, real_y, xlab = "X", ylab = "Y", type="l", ylim = c(-6,6))
points(x,y)

calc_natural_spline_basis <- function(x, knot, knot_max, dKm1){
  if(is.nan(knot)){
    return(x)
  }
  
  d1 = (relu((x-knot)^3) - relu(x-knot_max)^3)/(knot_max - knot)
  
  return(d1 - dKm1)
}

calc_natural_spline_mat <- function(x, knots){
  n = length(x)
  imat = matrix(rep(1,n), ncol = 1)
  imat = cbind(imat, x)
  
  knots = sort(knots)
  
  knot_max = knots[length(knots)]
  knot_second_max = knots[length(knots)-1]
  dKm1 = (relu((x-knot_second_max)^3) - relu(x-knot_max)^3)/(knot_max - knot_second_max)
  
  #browser()
  for(i in 1:(length(knots)-2)){
    imat = cbind(imat, calc_natural_spline_basis(x, knots[i], knot_max, dKm1))
  }
  
  return(imat)
}

natural_spline_regression <- function(x, y, nknots){
  knots = seq(min(x), max(x), length.out = (nknots+2))
  knots = knots[2:(length(knots)-1)]
  
  imat = calc_natural_spline_mat(x, knots)
  beta = ginv(crossprod(imat, imat)) %*% t(imat) %*% y
  y_hat = imat %*% beta
  
  return(list('knots'=knots, 'beta'=beta, 'y_hat'=y_hat))
}

ret = natural_spline_regression(x, y, nknots = 6)
lines(x, ret$y_hat, col="blue")
