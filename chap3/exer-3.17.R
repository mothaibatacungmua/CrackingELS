library(ElemStatLearn)

data("spam")

spam_dat = spam
num_col = ncol(spam_dat)
num_row = nrow(spam_dat)
spam_dat[,num_col] = as.character(levels(spam[,num_col]))[spam[,num_col]]
spam_dat[which(spam_dat[,num_col] == "spam"),num_col] = -1
spam_dat[which(spam_dat[,num_col] == "email"),num_col] = 1
spam_dat[,num_col] = as.numeric(spam_dat[,num_col])
names(spam_dat)[num_col] = "response"

# normal LS
library(boot)
set.seed(1000)
glm.fit = glm(response~0+., data = spam_dat)
cv.err = cv.glm(spam_dat, glm.fit, K=10)
print(cv.err$delta)

library(pls)
# PCR regression
pcr.fit = pcr(response~0+., ncomp = (num_col - 1), data = spam_dat, validation = "CV")
plot(RMSEP(pcr.fit), legendpos = "topright")
print(pcr.fit$coefficients)

# PLS regression
plsr.fit = plsr(response~0+., ncomp = (num_col - 1), data = spam_dat, validation = "CV")
plot(RMSEP(plsr.fit), legendpos = "topright")
print(pcr.fit$coefficients)


# Best subset
library(leaps)
library(MASS)
# Pls don't run it, if not you may wait until you die
bestsubset = regsubsets(response~0+., data=spam_dat, nvmax=(ncol - 1), really.big = TRUE)
summary(bestsubset)


library(glmnet)
#ridge regulariztion
X = model.matrix(response~0+., data = spam_dat)
y = spam_dat[,num_col]
ridge.cv = cv.glmnet(X, y)
plot(ridge.cv)
# lamda = exp(-6)
ridge.mod = glmnet(X, y, alpha = 0, lambda = exp(-6))
summary(ridge.mod$beta)

lasso.mod = glmnet(X, y, alpha = 1, lambda = exp(-6))
summary(lasso.mod$beta)


library(lars)
# draw lasso path
lars.fit = lars(X,y, type = "lar")
plot(lars.fit)
cv.lars(X,y,trace = TRUE, max.step = 200, plot.it = TRUE)
