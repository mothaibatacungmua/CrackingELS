library(ElemStatLearn)

attach(prostate)
#summary(prostate)
#head(prostate)

source('normalize_prostate_data.R')
dat = normalize_prostate_data(prostate)
df_train = dat[[1]]
df_test = dat[[2]]
ln = names(df_train)

# predictors matrix and response vector
X = model.matrix(~lcavol+lweight+age+lbph+svi+lcp+gleason+pgg45, data = df_train)
X = X[,-1]
y = df_train$lpsa

#using QR decomposition
QR <- qr(X)
Q <- qr.Q(QR)
R <- qr.Q(QR)

ln_w_i = ln[1:8]
active_set = c()
current_RSS = y
save_RSS = c(crossprod(current_RSS, current_RSS))
temp_Q = Q
n = ncol(temp_Q)
colnames(temp_Q) = ln_w_i

for (i in 1:n){
  pr = t(temp_Q) %*% current_RSS
  max_index = which.max(abs(pr))
  print(max_index)
  #add to the active set
  active_set = c(active_set, ln_w_i[max_index])
  
  #update current vector RSS
  current_RSS = current_RSS - crossprod(y, temp_Q[,max_index])*temp_Q[,max_index]
  
  #save value RSS
  save_RSS = c(save_RSS, crossprod(current_RSS, current_RSS))
  
  #remove the most related vector w.r.t current_RSS
  temp_Q = as.matrix(temp_Q[,-max_index])
  ln_w_i = ln_w_i[-max_index]
}

print(active_set)
print(save_RSS)

plot(0:n, save_RSS, 
     xlab = "Number of variables in the active set", 
     ylab = "Residual Sum-of-Squares",
     xlim=c(0,n), ylim=c(0,max(save_RSS)), 
     pch=21,bg="black")
lines(0:n, save_RSS)
adding = paste("+",active_set)
text(0:n, save_RSS,labels = c("NULL",adding), cex= 0.7,pos=1)

#using leaps package in short
library(leaps)
library(MASS)
forward_stepwise = regsubsets(lpsa~lcavol+lweight+age+lbph+svi+lcp+gleason+pgg45, data=df_train, nvmax=8, method="forward")
summary(forward_stepwise)

