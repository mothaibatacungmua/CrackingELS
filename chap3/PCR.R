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

# X don't have interpret
svd_pcr_regress <- function(X, y, maintain_var=0.01){
  svd_decomp = svd(X)
  D = svd_decomp$d
  V = svd_decomp$v
  
  k = 0
  for(i in 1:length(D)){
    if(sum(D[1:i])/sum(D) >= (1-maintain_var)){
      break
    }
    k = k + 1
  }
  
  Z = X %*% V[,1:k]
  beta = c(crossprod(Z,y)/diag(crossprod(Z,Z)))
  
  return(list(coefs=beta, Z=Z, V=V))
}

pcr_obj = svd_pcr_regress(X,y,0.1)
beta_pcr = pcr_obj$V %*% pcr_obj$coefs

# I want to use the NIPALS algorithm to solve PCR. It's the same as the SVD
# decomposition but it's the inside of PLS

library(base)

NIPALS <- function(X){
  save_X = X
  t_h = as.vector(X[,1])
  P = c()
  loop_count = 0
  
  while(loop_count < ncol(X)){
    while(TRUE){
      prev_t = t_h
      p_h = crossprod(save_X, t_h)/as.numeric(crossprod(t_h, t_h))
      p_h = p_h/sqrt(as.numeric(crossprod(p_h, p_h)))
      t_h = (X%*%p_h)
      
      if(sqrt(crossprod(prev_t-t_h, prev_t-t_h)) < 1e-8){
        save_X = save_X - t_h %*% t(p_h)
        P = c(P, p_h)
        loop_count = loop_count + 1
        break
      }
    }
  }
  
  P = matrix(P, nrow = ncol(X), byrow = FALSE)
  
  return(P)
}

P = NIPALS(X)
Z = X %*% P
beta = c(mean(y), crossprod(Z,y)/diag(crossprod(Z,Z)))
beta_pcr = P %*% beta

