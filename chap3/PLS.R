library(ElemStatLearn)

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

# The Partial Least Square is a extension of the NIPALS algorithm

PLS <- function(X, y){
  save_X = X
  u = y
  P = c()
  loop_count = 0
  T = c()
  
  while(loop_count < ncol(X)){
    w = crossprod(save_X,u)/norm(crossprod(save_X,u),type = "2")
    t = save_X %*% w
    p = crossprod(save_X, t)/as.numeric(crossprod(t, t))
    save_X = save_X - t%*%t(p)
    T = c(T, t)
    P = c(P, p)
    loop_count = loop_count + 1
  }
  
  P = matrix(P, nrow = ncol(X), byrow = FALSE)
  T = matrix(T, nrow = nrow(X), byrow = FALSE)
  
  return(list(P=P, T=T))
  #
  # With y is a matrix, we do the loop as the followings
  # u = Y[,j]
  # Loop:
  #   p = t(X)*u / ||t(X)*u||
  #   t = X*p
  #   q = t(Y)*t / ||t(Y)*t||
  #   u = Y*u
  # Until t stop changing
  #   X = X - t*t(p)
  #   Y = Y - u*t(q)
}

ret = PLS(X,y)
P = ret$P
T = ret$T

beta = crossprod(T,y)/diag(crossprod(Z,Z))
beta_pls = ginv(t(P)) %*% beta
