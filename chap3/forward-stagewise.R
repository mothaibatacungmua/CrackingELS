library(ElemStatLearn)

# Note: The nature of Forward Stagewise is the Steepest descent algorithms for L1-norm, 
# See more: chapter 9 in thet Convex Optimization book of Stephen Boyd

attach(prostate)
#summary(prostate)
#head(prostate)

source('normalize-prostate-data.R')
dat = normalize_prostate_data(prostate)
df_train = dat[[1]]
df_test = dat[[2]]
ln = names(df_train)

# predictors matrix and response vector
X = model.matrix(~lcavol+lweight+age+lbph+svi+lcp+gleason+pgg45, data = df_train)
X = X[,-1]
y = df_train$lpsa
colnames(X) = ln[1:8]

forward_stagewise <- function(X, y, epsilon, max_iter){
  p = ncol(X)  
  l1_norm = 0
  count_record = 0
  num_record = 30
  recorder = matrix(rep(0, p+1), ncol = (p+1))
  current_rss = y
  beta = rep(0, p)
  
  
  for(i in 1:max_iter){
    count_record = count_record + 1
    corr = crossprod(X, current_rss)
    prev_beta = beta
    argmax = which.max(abs(corr))
    
    #update beta
    beta[argmax] = beta[argmax] + epsilon*sign(corr[argmax])
    
    #update current_rss
    current_rss = current_rss - epsilon*sign(corr[argmax])*X[,argmax]
    
    l1_norm = sum(abs(beta))
    
    if(count_record == (max_iter/num_record) || prev_beta[argmax] == 0.0){
      recorder = rbind(recorder, c(beta, l1_norm))
      
      if(count_record == (max_iter/num_record)) count_record = 0
    }
    
  }  
  
  colnames(recorder) = c(colnames(X), "L1-norm")
  
  return(list(beta=beta, recorder=recorder))
}


fstaw_plot <- function(beta, recorder){
  p = length(beta)
  min_beta = min(recorder[,1:p])
  max_beta = max(recorder[,1:p])
  min_l1_norm = min(recorder[,p+1])
  max_l1_norm = max(recorder[,p+1])
  
  x_points = seq(min_l1_norm, max_l1_norm, length.out = nrow(recorder))
  y_points = rep(0, nrow(recorder))
  par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
  # create axis
  plot(x_points, y_points, type = "n", xlab = "L1-norm", ylab = "Coefficients",
       xlim = c(min_l1_norm, max_l1_norm), ylim = c(min_beta, max_beta))
  lines(x_points, y_points, type = "b", lwd=1.5, lty=2, pch="-")
  
  
  # plot coefficients
  colors = rainbow(p)
  l1_points = recorder[,p+1]
  for(i in 1:p){
    re_col = recorder[,i]
    
    index_not_null = 1
    while(re_col[index_not_null] == 0.0) index_not_null = index_not_null + 1
    
    points(l1_points[index_not_null-1], c(0), pch=21, bg=colors[i], cex=0.8)  
    
    lines(l1_points[(index_not_null-1):length(l1_points)], 
          re_col[(index_not_null-1):length(re_col)], 
          lwd=1.5, col=colors[i])
  }
  
  title("Forward Stagewise Regression")
  legend(max_l1_norm+0.5, 0.75, inset=c(-0.2,0), 
         legend=colnames(recorder[,1:p]), 
         cex = 0.5, col = colors, lty=1, title = "coefs")
  
}

epsilon = 0.001
max_iter = 30000
ret = forward_stagewise(X, y, epsilon, max_iter)
fstaw_plot(ret$beta, ret$recorder)
