# More details for this algorithm 
# See: http://www.stat.cmu.edu/~ryantibs/convexopt-S15/

soft_thresholding <- function(betas, lambda){
  betas[which(betas > lambda)] = betas[which(betas > lambda)] - lambda
  betas[which(betas < -lambda)] = betas[which(betas < -lambda)] + lambda
  
  betas[which(abs(betas) <= lambda)] = 0
  
  return(betas)
}

cost_fn <- function(X, y, betas, lambda){
  return((1/2)*crossprod(y - X %*% betas, y - X %*% betas) + lambda*sum(abs(betas[2:length(betas)])))  
}

prox_map <- function(betas, t, lambda){
  return(soft_thresholding(betas, t*lambda))
}

prox_lasso <- function(X, y, lambda, learning_rate = 0.001, max_iter=100000, thres=1e-8){
  monitor_cost = c()
  betas = rnorm(ncol(X),sd=5)/sum(sum(X))
  prev_value = 0
  #browser()
  next_value = cost_fn(X, y, betas, lambda)
  monitor_cost = c(monitor_cost, next_value)
  count_loop = 1
  while(count_loop < max_iter){
    #update beta
    betas = prox_map(betas + learning_rate*crossprod(X, y - X%*%betas), learning_rate, lambda)
    prev_value = next_value
    next_value = cost_fn(X, y, betas, lambda)
    
    monitor_cost = c(monitor_cost, next_value)
    
    if(abs(next_value - prev_value) < thres){
      break
    }
    #learning_rate = learning_rate/sqrt(count_loop)
    count_loop = count_loop + 1
    
    
  }
  
  return(list('beta'=betas,'monitor'=monitor_cost))
}

#Test case
set.seed(1)
X = cbind(floor(runif(5, 1, 5)), floor(runif(5, 1, 5)))
y = floor(runif(5, 1, 20))
X = cbind(rep(1,5),X) # adding interpret
ret = prox_lasso(X,y,0.6)
n_iter = length(ret$monitor)
plot(1:n_iter, ret$monitor, type = "l", main = "Proximal Gradient", xlab = "Iterations", ylab = "Cost")
