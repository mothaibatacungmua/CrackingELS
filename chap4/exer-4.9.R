ibrary(ElemStatLearn)

data("vowel.train")
data("vowel.test")
source("utils.R")
source("lda.R")
source("qda.R")

train_dat = vowel.train
test_dat = vowel.test

X_train = as.matrix(train_dat[,2:ncol(train_dat)])
y_train = as.vector(train_dat[,1])
X_test = as.matrix(test_dat[,2:ncol(test_dat)])
y_test = as.vector(test_dat[,1])

# tranforms y_train and y_test to the same as factors
y_train = convert_class_vector(y_train)
y_test = convert_class_vector(y_test)

lda_calc_error_wrapper(X_train, y_train, X_test, y_test)
qda_calc_error_wrapper(X_train, y_train, X_test, y_test)