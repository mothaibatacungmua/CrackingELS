library(ElemStatLearn)

attach(prostate)
#summary(prostate)
#head(prostate)

source('normalize-prostate-data.R')
dat = normalize_prostate_data(prostate)
df_train = dat[[1]]
df_test = dat[[2]]
ln = names(df_train)

# fit subset size k = 0
lm.fit = lm(predictors~+1, data=df_train)
rss = sum(lm.fit$residuals^2)
x_points = c(0) 
y_points = c(rss)
x_mins = c(0)
y_mins = c(rss)

n = 8
for (k in 1:n){
  list_combn = combn(n, k)
  x_mins = c(x_mins, k)
  y_mins = c(y_mins, 1000)
  
  for(i in 1:ncol(list_combn)){
    formu = paste(ln[list_combn[,i]],collapse = "+")
    formu = paste("lpsa~", formu)
    
    lm.fit = lm(formula = formu, data = df_train)
    rss = sum(lm.fit$residuals^2)
    x_points = c(x_points, k)
    y_points = c(y_points, rss)

    if(rss < y_mins[k+1]){
      y_mins[k+1] = rss
    }
  }
}

#hard-coded
y_points[1] = 100
plot(x_points, y_points, xlab = "Subset size k", ylab = "Residual Sum-of-Squares",
     ylim = c(0,100), xlim = c(0, n))

#hard-coded
y_mins[1] = 100
lines(x_mins, y_mins)

#using leaps package
library(leaps)
library(MASS)
bestsubset = regsubsets(lpsa~lcavol+lweight+age+lbph+svi+lcp+gleason+pgg45, data=df_train, nvmax=8)
summary(bestsubset)
