library(ElemStatLearn)
library(MASS)

# Note: The nature of Least Angle Regression a descent algorithms, 
# but its delta finding procedure and its line searching procedure
# are more efficient than the Forward Stagewise

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

find_min_positive <- function(v1, v2){
  vl = 0
  
  if(v1 > 0) vl = v1
  if(v2 > 0) vl = v2
  if(v1 > 0 && v2 >0) vl = min(v1, v2)
  
  return(vl)
}

look_ahead <- function(c_max, vec_a_k, inactive_set, corr){
  # The following code due to the solution of the exercie 3.25(ESL), it's the look-ahead properties of LAR
  
  if(length(inactive_set) == 0)
    return(list(min_line=NULL, min_value=1))
  
  #print(inactive_set)
  m_i = inactive_set[1]
  v1 = (c_max - corr[m_i])/(c_max - vec_a_k[m_i])
  v2 = (c_max + corr[m_i])/(c_max + vec_a_k[m_i])
  vl = find_min_positive(v1, v2)
  
  if(length(inactive_set) >= 2){ 
    #fuck you, R script. The looping also occur when even length(inactive_set) < 2
    for(i in 2:length(inactive_set)){
      index = inactive_set[i]
      v1 = (c_max - corr[index])/(c_max - vec_a_k[index])
      v2 = (c_max + corr[index])/(c_max + vec_a_k[index])
      min_vl = find_min_positive(v1, v2)
      
      if(min_vl < vl){
        m_i = index
        vl = min_vl
      }
    } 
  }
  
  return(list(min_line=m_i, min_value=vl))
} 

find_drop_index <- function(current_beta, active_set, direction, vl){
  steps = -(current_beta[active_set]/direction)
  drop_var = -1
  steps[which(steps <= 0)] = 2e308
  min_hit_zero_index = which.min(steps)
  
  if(steps[min_hit_zero_index] <= vl){
    drop_var = active_set[min_hit_zero_index]
    vl = steps[min_hit_zero_index]
  }
  
  return(list(drop_var=drop_var, vl=vl))
}

lars_lasso <- function(X, y){
  # init input
  ln_w_i = colnames(X)
  p = ncol(X)
  current_rss = y
  inactive_set = 1:p
  active_set = c()
  recorder = matrix(rep(0, p+1), nrow = 1)
  
  # step 1: find the most correlated variable with r
  corr = t(X) %*% current_rss
  argmax = which.max(abs(corr))
  active_set = c(argmax)
  inactive_set = inactive_set[-which(inactive_set==argmax)]
  beta = rep(0, p)
  last_update = FALSE
  
  # step 2: moving beta and adding the remain variables
  while(length(active_set) != p || !last_update){
    active_X = X[,active_set]
    c_max = abs(corr[argmax])
    #browser()
    # caculate delta
    direction = ginv(crossprod(active_X, active_X)) %*% t(active_X) %*% current_rss
    fit_direction = active_X %*% direction
    
    # line search
    vec_a_k = t(X) %*% fit_direction
    ret = look_ahead(c_max, vec_a_k, inactive_set, corr)
    m_i = ret$min_line
    vl = ret$min_value
    
    last_update = is.null(m_i)
    
    # lars modified for finding the lasso path
    f = find_drop_index(beta, active_set, direction, vl)
    drop_var = f$drop_var
    vl = f$vl
    
    # update beta and track L1-norm
    if(drop_var != -1){
      inactive_set = c(inactive_set, drop_var)
      drop_index = which(active_set == drop_var)
      active_set = active_set[-drop_index]
      beta[drop_var] = 0
      beta[active_set] = beta[active_set] + vl*direction[-drop_index]
    }else{
      beta[active_set] = beta[active_set] + vl*direction
      # update active set
      active_set = c(active_set, m_i)
      inactive_set = inactive_set[-which(inactive_set==m_i)]
    }
    
    l1_norm = sum(abs(beta))
    recorder = rbind(recorder, c(beta, l1_norm))
    
    # update current rss
    current_rss = current_rss - vl*fit_direction
    
    
    # update correlation vector
    corr = t(X) %*% current_rss
    
  } 
  
  colnames(recorder) = c(ln_w_i, "L1-norm")
  
  return(list(beta=beta, recorder=recorder))
}

lasso_lars_plot <- function(beta, recorder){
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
  
  title("Lasso Path")
  
  legend(max_l1_norm*1.2, max_beta, inset=c(-0.2,0), 
         legend=colnames(recorder[,1:p]), 
         cex = 0.5, col = colors, lty=1, title = "coefs")
}

ret = lars_lasso(X,y)
print(ret$beta)
lasso_lars_plot(ret$beta, ret$recorder)
