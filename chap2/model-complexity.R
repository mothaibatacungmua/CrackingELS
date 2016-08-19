#based on code on the statistic course: http://www.stat.rice.edu/~gallen/stat640.html
library(class)

set.seed(1)
n = 300
sig = 1.5
mu1 = c(7,2); mu2 = c(2,7); mu3 = c(7,7)

x1 = t(matrix(mu1, 2, n/3)) + matrix(rnorm(n)*sig, n/3, 2)
x2 = t(matrix(mu2, 2, n/3)) + matrix(rnorm(n)*sig, n/3, 2)
x3 = t(matrix(mu3, 2, n/3)) + matrix(rnorm(n)*sig, n/3, 2)
x.train = rbind(x1, x2, x3)

set.seed(100)
x1 = t(matrix(mu1, 2, n/3)) + matrix(rnorm(n)*sig, n/3, 2)
x2 = t(matrix(mu2, 2, n/3)) + matrix(rnorm(n)*sig, n/3, 2)
x3 = t(matrix(mu3, 2, n/3)) + matrix(rnorm(n)*sig, n/3, 2)
x.test = rbind(x1, x2, x3)

y.train = c(rep(1, n/3), rep(2, n/3), rep(3, n/3))
Y = factor(y.train)

loop = 100
train_error = rep(0, loop)
test_error = rep(0, loop)

for (i in 1:loop){
  knn.pred = knn(train = x.train, test = x.train, cl = y.train, k = i)
  train_error[i] = mean(knn.pred == y.train)
  knn.pred = knn(train = x.train, test = x.test, cl = y.train, k = i)
  test_error[i] = mean(knn.pred == y.train)
}

k_s = 1:loop
library(ggplot2)
df <- data.frame(k_s,train_error,test_error)
ggplot(df, aes(k_s)) + geom_line(aes(y=train_error), colour = "red") + geom_line(aes(y=test_error), colour = "green")
