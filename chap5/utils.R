gen_data <- function(){
  x = seq(from=0, to=1, length.out = 50)
  n = length(x)
  y = 4*sin(-6*x+2) + rnorm(n, 0, 2)
  real_y = 4*sin(-6*x+2)
  return(list('x'=x,'y'=y, 'real_y'=real_y))
}

relu <- function(x){
  z = rep(0, length(x))
  ind = which(x > z)
  z[ind] = x[ind]
  
  return(z)
}