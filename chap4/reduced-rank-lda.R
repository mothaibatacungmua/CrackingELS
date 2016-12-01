library(ElemStatLearn)

data("vowel.train")
data("vowel.test")
source("utils.R")
source("lda.R")

train_dat = vowel.train
test_dat = vowel.test

X_train = as.matrix(train_dat[,2:ncol(train_dat)])
y_train = as.vector(train_dat[,1])
X_test = as.matrix(test_dat[,2:ncol(test_dat)])
y_test = as.vector(test_dat[,1])

# tranforms y_train and y_test to the same as factors
y_train = convert_class_vector(y_train)
y_test = convert_class_vector(y_test)

W = calc_within_class(X_train, y_train)
B = calc_between_class(X_train, y_train)
calc_discriminant_coordinates <- function(B, W){
  svd_W = svd(W) 
  squared_W = svd_W$u %*% diag(svd_W$d^(-1/2)) %*% t(svd_W$v)
  Vl = svd(t(squared_W) %*% B %*% squared_W)$u
  
  return(squared_W %*% Vl)
}

calc_discriminant_variables <- function(X, discr_coords){
  return(X %*% discr_coords)
}

# At here, I make a mistake, one should sphere data before reducing because
# the covariance matrix of variables in X-sphered -> X-sphered-reduced not change
# and equals I.
discr_coords = calc_discriminant_coordinates(B, W)
X_trans = calc_discriminant_variables(X_train, discr_coords)
X_test_trans = calc_discriminant_variables(X_test, discr_coords)

# plot results
nl = nlevels(as.factor(y_train))
cols = rainbow(nl)
ycols = y_train
for(i in 1:length(cols)){
  ycols[which(y_train==i)] = cols[i]
}
plot(X_trans[,1], X_trans[,2], col=ycols,bg=ycols, pch=21)


lda_calc_error_wrapper(X_trans[,1:2] ,y_train, X_test_trans[,1:2], y_test)
