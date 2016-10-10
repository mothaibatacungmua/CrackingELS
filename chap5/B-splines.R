library(MASS)

source('./utils.R')

# the knots parameter is included two boundary knots.
calc_b_spline <- function(x, knots, max_order, order, idx){
  #browser()
  internal_func <- function(x, aug_knots, m, i){
    #browser()
    B = rep(0, length(x))
    if(m == 1){
      if(aug_knots[i+1] == max(aug_knots)){
        indx = which(x >= aug_knots[i] & x <= aug_knots[i+1])  
      }else{
        indx = which(x >= aug_knots[i] & x < aug_knots[i+1])
      }
      
      B[indx] = 1
      
      return(B)
    }
    
    fzero = TRUE
    for(j in i:(i+m)){
      if(aug_knots[j]!=aug_knots[i]){
        fzero = FALSE
        break
      } 
    }
    
    if(fzero) {
      return(B)
    }
    
    B_left = internal_func(x, aug_knots, m-1, i)
    B_right = internal_func(x, aug_knots, m-1, i+1)
    
    if(aug_knots[i] != aug_knots[i+m-1]){
      B = B + (x - aug_knots[i])/(aug_knots[i+m-1] - aug_knots[i]) * B_left
    }
    
    if(aug_knots[i+1] != aug_knots[i+m]){
      B = B + (aug_knots[i+m] - x)/(aug_knots[i+m] - aug_knots[i+1]) * B_right
    }
    
    return(B)
  }
  
  bmax = max(knots)
  bmin = min(knots)
  aug_left = rep(bmin, max_order-1)
  aug_right = rep(bmax, max_order-1)
  aug_knots = c(aug_left, knots, aug_right)
  
  return(internal_func(as.vector(x), aug_knots, order, idx))
}

b_spline <- function(x, knots, max_m, m, i){
  #browser()
  if(max(knots) >= max(x) | min(knots) <= min(x)){
    return(NULL)
  } 
  
  equal = FALSE
  for(k in 1:length(knots)-1){
    findk = which(knots == knots[k])
    if(length(findk) >= 2){
      equal = TRUE
    } 
  }
  
  if(equal){
    return(NULL)
  }
  
  added_boundary_knots = c(min(x), knots, max(x))
  y = calc_b_spline(x, added_boundary_knots, max_m, m, i)
  
  return(y)
}

dat = gen_data()
x = dat$x
y = dat$y
real_y = dat$real_y

par(mfrow = c(2, 1))
q = quantile(x)
colors = rainbow(7)
for(i in 1:7){
  yi = b_spline(x, q[2:(length(q)-1)], 4, 4, i)
  if(i == 1){
    plot(x, yi, type="l", col=colors[i], main = "B-Splines", xlab="x", ylab="y")
    next
  }
  lines(x, yi, type="l", col=colors[i])
}


b_spline_regression <- function(x, y, knots, order){
  num_basis = length(knots) + order
  N = length(x)
  b_mat = matrix(rep(0, N*num_basis), nrow = N)
  
  for(i in 1:num_basis){
    b_mat[,i] = b_spline(x, knots, order, order, i)
  }
  
  beta = ginv(crossprod(b_mat, b_mat)) %*% t(b_mat) %*% y
  y_hat = b_mat %*% beta
  
  return(list('beta'=beta, 'y_hat'=y_hat))
}

plot(x, real_y, xlab = "X", ylab = "Y", type="l", ylim = c(-6,6), main = "B-Splines Regression")
points(x,y)

ret = b_spline_regression(x, y, knots = q[2:(length(q)-1)], order=4)
lines(x, ret$y_hat, col="red")