# barrier method for solving the hyperplane optimization problem
# LP problem:
#   min_{beta, beta0} 1/2||beta||^2
#   subject to yi(t(xi)*beta + beta0) >= 1, i=1,...,N

library(MASS)
library(base)
source('utils.R')
source('gen-data.R')

dat = gen_two_class_data()
X = dat$X
y = dat$y

X = X[,2:3]
plot(X[,1], X[,2], col=as.numeric(y)+2, xlab = "X1", ylab = "X2", bg=as.numeric(y)+2, pch=21, main = "Separating Hyperplane")


phase_I_cost <- function(X, y, t, s, beta, beta0){
  f = 1 - y*(X %*% beta + beta0)-s
  greater_zero = which(f>0)
  
  if(length(greater_zero) > 0){
    return(1/0)
  }
  
  return(t*s + sum(-log(-f)))
}

phase_I_grad <- function(X, y, t, s, beta, beta0){
  f = 1 - y*(X %*% beta + beta0)-s
  grad_s = t + sum(1/f)
  grad_beta0 = sum(y/f)
  grad_beta = t(X*y) %*% (1/f)
  
  return(list('s'=grad_s,'beta0'=grad_beta0,'beta'=grad_beta))
}

phase_I_hessian <- function(X, y, t, s, beta, beta0){
  f = 1 - y*(X %*% beta + beta0)-s
  hess_s = sum(1/f^2)
  hess_beta0 = hess_s
  hess_beta = crossprod(
    X/as.vector(f), X/as.vector(f)
  )
  
  return(list('s'=hess_s,'beta0'=hess_beta0,'beta'=hess_beta))
}

phase_I_line_search <- function(X, y, t, s, beta, beta0){
  alpha = 0.25
  gamma = 0.7
  l = 1
  grad = phase_I_grad(X, y, t, s, beta, beta0)
  hessian = phase_I_hessian(X, y, t, s, beta, beta0)
  
  ds = -as.numeric((1/hessian$s)*grad$s)
  dbeta = -ginv(hessian$beta) %*% as.vector(grad$beta)
  dbeta0 = -as.numeric((1/hessian$beta0) * grad$beta0)
  
  #using backtracking
  while(TRUE){
    cost =  phase_I_cost(X, y, t, s+l*ds, beta+l*dbeta, beta0+l*dbeta0)
    if(is.infinite(cost)){
      while(is.finite(cost)){
        l = gamma*l
        cost =  phase_I_cost(X, y, t, s+l*ds, beta+l*dbeta, beta0+l*dbeta0)
      }
      
      return(l)
    }
    
    bound = phase_I_cost(X, y, t, s, beta, beta0) + alpha*l*(grad$s*ds + t(grad$beta)%*%dbeta + grad$beta0*beta0)
    if(as.numeric(cost) < as.numeric(bound)){
      break
    }
    
    l = gamma*l
  }
  
  return(l)
}

# Phase I, find a feasible solution. 
phase_I_find_t_op <- function(X, y, t, s, beta, beta0, max_iter=40){
    count = 0
    #browser()
    while(count < max_iter){
      grad = phase_I_grad(X, y, t, s, beta, beta0)
      hessian = phase_I_hessian(X, y, t, s, beta, beta0)
      
      l = phase_I_line_search(X, y, t, s, beta, beta0)
      #update s, beta, beta0
      s = s - l*as.numeric(1/hessian$s * grad$s)
      beta0 = beta0 - l*as.numeric(1/hessian$beta0 * grad$beta0)
      beta = beta - l*ginv(hessian$beta) %*% as.vector(grad$beta)
      count = count + 1
      
      if(s < 0){
        break
      } 
    }
    
    return(list('s'=s,'beta0'=beta0,'beta'=beta))
}

phase_I <- function(X, y){
  nl = ncol(X)
  beta = rnorm(nl,sd=3)
  beta0 = 100
  
  f = 1 - y*(X %*% beta + beta0)
  
  # add one at here to avoid dividing zero
  s = max(f)+1
  
  if(s <= 0){
    return(list('beta':beta, 'beta0':beta0))
  }
  
  epsilon = 0.00000001
  m= 32
  t = 1
  u = 3

  while(m/t > epsilon){
    min_t = phase_I_find_t_op(X, y, t, s, beta, beta0)
    s = min_t$s
    beta0 = min_t$beta0
    beta = min_t$beta
    if(s < 0){
      break
    }
    
    t = t*u
  }
  
  return(list('beta'=beta, 'beta0'=beta0))
}

#using backtracking
phase_II_cost <- function(X, y, t, beta, beta0){
  f = 1 - y*(X %*% beta + beta0)
  greater_zero = which(f>=0)
  
  if(length(greater_zero) > 0){
    return(1/0)
  }
  
  return(1/2*sum(beta^2) + sum(-log(-f)))
}

phase_II_grad <- function(X, y, t, beta, beta0){
  f = 1 - y*(X %*% beta + beta0)
  
  grad_beta0 = sum(y/f)
  grad_beta = t*beta + t(X*y) %*% (1/f)
  
  return(list('beta0'=grad_beta0,'beta'=grad_beta))
}

phase_II_hessian <- function(X, y, t, beta, beta0){
  f = 1 - y*(X %*% beta + beta0)
  
  hess_beta0 = sum(1/f^2)
  hess_beta = t*diag(ncol(X)) + crossprod(
    X/as.vector(f), X/as.vector(f)
  )
  
  return(list('beta0'=hess_beta0,'beta'=hess_beta))
}

phase_II_line_search <- function(X, y, t, beta, beta0){
  alpha = 0.5
  gamma = 0.7
  l = 1
  grad = phase_II_grad(X, y, t, beta, beta0)
  hessian = phase_II_hessian(X, y, t, beta, beta0)
  
  dbeta = -ginv(hessian$beta) %*% as.vector(grad$beta)
  dbeta0 = -as.numeric((1/hessian$beta0)*grad$beta0)
  

  cost =  phase_II_cost(X, y, t, beta+l*dbeta, beta0+l*dbeta0)
  if(is.infinite(cost)){
    while(is.infinite(cost)){
      l = gamma*l
      cost =  phase_II_cost(X, y, t, beta+l*dbeta, beta0+l*dbeta0)
    }
    
    #bound = phase_II_cost(X, y, t, beta, beta0) + alpha*l*(t(grad$beta)%*%dbeta + grad$beta0*beta0)
    #if(as.numeric(cost) <= as.numeric(bound)){
    #  break
    #}
    
    #l = gamma*l
  }
    
  return(l)
}

phase_II_find_t_op <- function(X, y, t, beta, beta0, max_iter=50){
  count = 0
  #browser()
  while(count < max_iter){
    f = 1 - y*(X %*% beta + beta0)
    
    grad = phase_II_grad(X, y, t, beta, beta0)
    hessian = phase_II_hessian(X, y, t, beta, beta0)
    l = phase_II_line_search(X, y, t, beta, beta0)
    
    #update s, beta, beta0
    beta0 = beta0 - l*as.numeric(1/hessian$beta0 * grad$beta0)
    beta = beta - l*ginv(hessian$beta) %*% as.vector(grad$beta)
    count = count + 1
  }
  
  return(list('beta0'=beta0,'beta'=beta))
}

phase_II <- function(X, y, beta, beta0){
  epsilon = 0.001
  m= 1000
  t = 1.5
  u = 1.001
  
  while(m/t > epsilon){
    min_t = phase_II_find_t_op(X, y, t, beta, beta0)
    beta0 = min_t$beta0
    beta = min_t$beta
    t = t*u
  }
  
  return(list('beta'=beta, 'beta0'=beta0))
}

min_t = phase_I(X,y)
f= 1 - y*(X %*% min_t$beta + min_t$beta0)
nonzero = which(f>0)
if(length(nonzero) > 0){
  print('There arent exist a solution')
}else{
  result = phase_II(X,y,min_t$beta, min_t$beta0)
  beta = result$beta/sqrt(sum(result$beta^2))
  beta0 = result$beta0/sqrt(sum(result$beta^2))
  abline(-beta0/beta[2], -beta[1]/beta[2])
}