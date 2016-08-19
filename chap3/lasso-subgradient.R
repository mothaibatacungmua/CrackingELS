# More details for this algorithm 
# See: http://www.stat.cmu.edu/~ryantibs/convexopt-S15/

calc_subgradient <- function(X, y, betas, lambda){
  subgrad = rep(0, length(betas))
  subgrad[1] = -crossprod(X[,1],(y-X%*%betas))
  
  for(i in 2:length(betas)){
    if(betas[i] != 0){
      subgrad[i] = -crossprod(X[,i],(y-X%*%betas)) + lambda*sign(betas[i])
      next
    } 
    
    subgrad[i] = -crossprod(X[,i],(y-X%*%betas)) - lambda
  }
  
  return(subgrad)
}

cost_fn <- function(X, y, betas, lambda){
  return((1/2)*crossprod(y - X %*% betas, y - X %*% betas) + lambda*sum(abs(betas[2:length(betas)])))  
}

subgrad_lasso <- function(X, y, lambda, learning_rate = 0.001, max_iter=100000){
 betas = rnorm(ncol(X))
 keep_min = betas
 
 monitor_cost = c()
 
 prev_value = 0
 next_value = cost_fn(X, y, betas, lambda)
 count_loop = 1
 step = learning_rate
 monitor_cost = c(monitor_cost, next_value)
 
 while(count_loop < max_iter){
   #update beta
   betas = betas - step*calc_subgradient(X, y, betas, lambda)
   prev_value = next_value
   v = cost_fn(X, y, keep_min, lambda)
   next_value = cost_fn(X, y, betas, lambda)
   
   monitor_cost = c(monitor_cost, next_value)
   
   if(next_value < v){
     keep_min = betas
   }
   count_loop = count_loop + 1
   
 }
 
 return(list('beta'=keep_min,'monitor'=monitor_cost))
}

#Test case
set.seed(1)
X = cbind(floor(runif(5, 1, 5)), floor(runif(5, 1, 5)))
y = floor(runif(5, 1, 20))
X = cbind(rep(1,5),X) # adding interpret
ret = subgrad_lasso(X,y,0.6)
plot(1:10000, ret$monitor[1:10000], type = "l", main = "Subgradient", xlab = "Iterations", ylab = "Cost")

