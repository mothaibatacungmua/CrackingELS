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
df = data.frame(cbind(predictors, response))
ln = names(df)

# fit subset size k = 0
lm.fit = lm(predictors~+1, data=df)
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
    
    lm.fit = lm(formula = formu, data = df)
    rss = sum(lm.fit$residuals^2)
    x_points = c(x_points, k)
    y_points = c(y_points, rss)

    if(rss < y_mins[k+1]){
      print(rss)
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
