library(ElemStatLearn)
library(ggplot2)
library(reshape2)
library(MASS)
library(mvtnorm)

gaussian_func <- function (x, mean = 0, sd = 1) {
  return(dnorm(x, mean, sd))
}

mixture_gaussian <- function (x, mean_vec, sd_vec, coefs){
  ret = 0
  #normalize coefs
  coefs = coefs/sqrt(t(coefs)*coefs)
  for(i in 1:length(coefs)){
    ret = ret + coefs[i]*gaussian_func(x, mean_vec[i], sd_vec[i])
  }
  return(ret)
}


set.seed(1)
#
# Draw gaussian
#
x = seq(-6,6,length=200)
y1 = gaussian_func(x, -2, 1)
y2 = gaussian_func(x, 2, 1)
df <- data.frame(x, y1, y2)
df <- melt(df, id="x")
ggplot(df, aes(x = x, y = value, colour = variable)) + geom_line()

#
# Draw mixture gaussian
#
y = mixture_gaussian(x, c(-2, 2), c(1, 1), c(0.5, 0.5))
df <- data.frame(x, y)

#
# Generate multivarite gaussian functions
#
set.seed(1)
gen_random_sigma <- function (x){
  sigma = matrix(rep(0,4), ncol = 2)
  sigma[1,1] = runif(1, min = 0, max = 1)
  sigma[2,2] = runif(1, min = 0, max = 1)
  
  return(sigma)
}
v_sigma = matrix(sapply(1:10, gen_random_sigma), ncol = 2, byrow = T)
v_mean_0 = matrix(runif(20, min = -0.5, max = 0), ncol = 2)
v_mean_1 = matrix(runif(20, min = 0, max = 0.3), ncol = 2)

# function is used to generate random samples for each class
gen_sample_class <- function(n, vec_mu){
    ret = matrix(rep(0, 2*n), ncol = 2)
    for(i in 1:n){
      s = c(0, 0)
      for(k in 1:10){
        s = s + mvrnorm(1, mu = vec_mu[k,], Sigma = v_sigma[(k*2-1):(k*2),])
      }
      ret[i,] = s
    }
    
    return(ret)
}

# function is used to caculate log(density) function for a predictor with respect to a class
d_mix <- function (x, vec_mu){
  ret = 0
  for(k in 1:10){
    ret = ret + dmvnorm(x, mean = vec_mu[k,], sigma = v_sigma[(k*2-1):(k*2),])
  }
  
  return(log(ret))
}

N = 200
n0 = 100
n1 = 100

#
# Generate for class 0
#
y0 = rep(0, n0)
X0 = gen_sample_class(n0, v_mean_0)

#
# Generate for class 1
#
y1 = rep(1, n1)
X1 = gen_sample_class(n1, v_mean_1)

y = c(y0, y1)
X = rbind(X0, X1)

plot(X[,1], X[,2], col = as.numeric(y)+1, xlab = "X1", ylab = "x2")

gx1 = seq(min(X[,1]), max(X[,1]), length = 200)
gx2 = seq(min(X[,2]), max(X[,2]), length = 200)
pos = expand.grid(gx1, gx2)
d_ratios = apply(pos, 1, function(x){
  return(d_mix(x, v_mean_0) - d_mix(x, v_mean_1))
})

pr <- matrix(d_ratios, length(gx1), length(gx2))
contour(gx1, gx2, pr, levels = seq(-0.01,0.01,length = 10), axes = F, labels="", xlab="", ylab="", main=
          "Bayesian desicion boundary")
points(X[,1], X[,2], col = as.numeric(y)+1, xlab = "X1", ylab = "x2")
points(pos, pch=".", cex=1.2, col=ifelse(d_ratios < 0, "coral", "cornflowerblue"))
box()
