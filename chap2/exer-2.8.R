library(ElemStatLearn)
library(ggplot2)
library(class)

data("zip.train")
data("zip.test")

train.d = zip.train[which(zip.train[,1] == 2 | zip.train[,1] == 3), ]
test.d = zip.test[which(zip.test[,1] == 2 | zip.test[,1] == 3), ]

train.d[which(train.d[,1] == 2), 1] = -1
train.d[which(train.d[,1] == 3), 1] = 1

test.d[which(test.d[,1] == 2), 1] = -1
test.d[which(test.d[,1] == 3), 1] = 1

knn_train_err = c()
knn_test_err = c()

for(i in c(1,3,5,7,15)){
  knn.pred = knn(train.d[,2:ncol(train.d)], train.d[,2:ncol(train.d)], train.d[,1], k = i)
  knn_train_err = c(knn_train_err, 1 - mean(knn.pred == train.d[,1]))
  
  knn.pred = knn(train.d[,2:ncol(train.d)], test.d[,2:ncol(test.d)], train.d[,1], k = i)
  knn_test_err = c(knn_test_err, 1 - mean(knn.pred == test.d[,1]))
}

k_s = c(1,3,5,7,15)
df <- data.frame(k_s, knn_train_err, knn_test_err)
ggplot(df, aes(k_s)) + 
  geom_line(aes(y=knn_train_err, color = "Train error")) + 
  geom_line(aes(y=knn_test_err, color = "Test error")) +
  scale_colour_manual(values = c("Train error"="red", "Test error"="green"))+
  labs(title = "K-Nearest-Neighbors")


lm.fit = lm(train.d[,1]~0+train.d[,2:ncol(train.d)])
coefs = coef(lm.fit)
lm.pred = test.d[,2:ncol(test.d)] %*% coefs
lm.pred[which(lm.pred <= 0)] = -1
lm.pred[which(lm.pred > 0)] = 1

lm.error = 1 - mean(lm.pred == test.d[,1])

