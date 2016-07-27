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

#using QR decomposition
QR <- qr(X)
Q <- qr.Q(QR)
R <- qr.Q(QR)

ln_w_i = c("Interpret", ln[1:8])
removed_set = c()
current_RSS = y - Q %*% crossprod(Q, y)
save_RSS = c(crossprod(current_RSS, current_RSS))
temp_Q = Q
p = ncol(temp_Q)
N = nrow(temp_Q)
colnames(temp_Q) = ln_w_i

for(i in 1:p){
  #compute squared z-score
  z_score = c()
  rss_1 = y - temp_Q %*% crossprod(temp_Q, y)
  for(j in 1:ncol(temp_Q)){
    t = as.matrix(temp_Q[,-j])
    rss_0 = y - t %*% crossprod(t, y)
    #dont need to multiply with (ncol(temp_Q) - ncol(t))/(N - ncol(temp_Q) - 1)
    score = (crossprod(rss_0, rss_0) - crossprod(rss_1, rss_1))/crossprod(rss_1, rss_1)
    z_score = c(z_score, score)
  }
  
  argmax = which.min(z_score)
  removed_set = c(removed_set, ln_w_i[argmax])
  temp_Q = as.matrix(temp_Q[,-argmax])
  ln_w_i = ln_w_i[-argmax]
  current_RSS = y - temp_Q %*% crossprod(temp_Q, y)
  save_RSS = c(save_RSS, crossprod(current_RSS, current_RSS))
}

print(removed_set)
print(save_RSS)

plot(0:9, save_RSS, 
     xlab = "Number of variables in the removed set", 
     ylab = "Residual Sum-of-Squares",
     xlim=c(0,9), ylim=c(0,max(save_RSS)), 
     pch=21,bg="black")
lines(0:9, save_RSS)
adding = paste("-",removed_set)
text(0:9, save_RSS,labels = c("NULL", adding), cex= 0.7,pos=1)

#using leaps package in short
library(leaps)
library(MASS)
backward_stepwise = regsubsets(lpsa~lcavol+lweight+age+lbph+svi+lcp+gleason+pgg45, data=df, nvmax=8, method="backward")
summary(backward_stepwise)

#it can have a little different with my code due to the normalizing