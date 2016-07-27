# Demostrate caculate beta_hat by using gram-schmidt procedure
library(MASS)

attach(Boston)

lm.fit = lm(medv ~ crim + age + dis + lstat, data = Boston)
summary(lm.fit)

coefs <- lm.fit$coefficients
coefs


# gram-schmidt procedure phase 1, compute z_1, z_2, ..., z_{p+1}
N = nrow(Boston)
p = 4+1
X = cbind(rep(1, N), Boston$crim, Boston$age, Boston$dis, Boston$lstat)
y = Boston$medv

gram_schmidt <- function (X, y){
  Z = matrix(rep(0, N*p), ncol = p)
  Z[,1] = X[,1]
  for (j in 2:p){
    Z[,j] = X[,j]
    print(j)
    for (l in 1:(j-1)){
      cof = (t(X[,j]) %*% Z[,l]) / (t(Z[,l]) %*% Z[,l])
      Z[,j]  = Z[,j] - Z[,l] * cof
    }
  }
  
  return(Z)
}

# reverse to compute beta_hat
calc_beta <- function(X, y, Z){
  beta_hat = rep(0, p)
  for (j in p:1){
    accum = y
    if(j < p){
      for (k in (j+1):p){
        accum = accum - beta_hat[k]*X[,k]
      }
    }
    
    beta_hat[j] = (t(accum) %*% Z[,j]) / (t(Z[,j]) %*% Z[,j])
  }
  
  return(beta_hat)
}

Z = gram_schmidt(X, y)
beta_hat = calc_beta(X, y, Z)

(beta_hat)
(coefs)

QR <- qr(X)
Q <- qr.Q(QR)
R <- qr.R(QR)
backsolve(R, crossprod(Q,y))
