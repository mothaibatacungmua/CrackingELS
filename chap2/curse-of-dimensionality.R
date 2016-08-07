library(class)
set.seed(1)
n = 300
p0 = 1000
sig = 1.5
mu1 = c(2,7); mu2 = c(7,2); mu3 = c(7,7)
x1 = t(matrix(mu1, 2, n/3)) + matrix(rnorm(n)*sig, n/3, 2)
x2 = t(matrix(mu2, 2, n/3)) + matrix(rnorm(n)*sig, n/3, 2)
x3 = t(matrix(mu3, 2, n/3)) + matrix(rnorm(n)*sig, n/3, 2)
x.train = cbind(rbind(x1, x2, x3), matrix(rnorm(n*p0),n,p0))

x1 = t(matrix(mu1, 2, n/3)) + matrix(rnorm(n)*sig, n/3, 2)
x2 = t(matrix(mu2, 2, n/3)) + matrix(rnorm(n)*sig, n/3, 2)
x3 = t(matrix(mu3, 2, n/3)) + matrix(rnorm(n)*sig, n/3, 2)
x.test = cbind(rbind(x1, x2, x3), matrix(rnorm(n*p0),n,p0))

y.train = c(rep(1, n/3), rep(2, n/3), rep(3, n/3))
Y = factor(y.train)

knn.pred = knn(train = x.train, test = x.test, cl = y.train, k = 5)
mean(knn.pred == y.train)
