library(ElemStatLearn)

attach(prostate)
summary(prostate)
head(prostate)

training_data = prostate[train, ]
test_data = prostate[train == F, ]

# normalize predictors
predictors = training_data[,1:8]
predictors = scale(predictors)

attr(predictors, 'scaled:center')
attr(predictors, 'scaled:scale')

# extract response and create data frame
lpsa = training_data[,9]
df = data.frame(cbind(predictors, lpsa))
ln = names(df)

# predictors matrix and response vector
X = model.matrix(~lcavol+lweight+age+lbph+svi+lcp+gleason+pgg45, data = df)
y = df$lpsa

#using QR decomposition
QR <- qr(X)
Q <- qr.Q(QR)
R <- qr.Q(QR)

ln_w_i = c("Interpret", ln[1:8])
active_set = c()
current_RSS = y
save_RSS = c(crossprod(current_RSS, current_RSS))
temp_Q = Q
n = ncol(temp_Q)
colnames(temp_Q) = ln_w_i

for (i in 1:n){
  pr = t(temp_Q) %*% current_RSS
  max_index = which.max(pr)
  
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

plot(0:9, save_RSS, 
     xlab = "Number of variables in the active set", 
     ylab = "Residual Sum-of-Squares",
     xlim=c(0,9), ylim=c(0,max(save_RSS)), 
     pch=21,bg="black")
lines(0:9, save_RSS)
adding = paste("+",active_set)
text(0:9, save_RSS,labels = c("NULL", adding), cex= 0.7,pos=1)

#using leaps package in short
library(leaps)
library(MASS)
forward_stepwise = regsubsets(lpsa~lcavol+lweight+age+lbph+svi+lcp+gleason+pgg45, data=df, nvmax=8, method="forward")
summary(forward_stepwise)
