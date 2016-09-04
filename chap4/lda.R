library(MASS)

source('gen-data.R')
source('utils.R')

train_data = gen_train_data()
test_data = gen_test_data()
X = train_data$X
y = train_data$y 

plot(X[,1], X[,2], col=as.numeric(y)+1, xlab = "X1", ylab = "X2", bg=as.numeric(y)+1, pch=21)

# using normal computing
lda_linear_discr <- function(x, s, mu, pi_c){
  return(t(x) %*% ginv(s) %*% mu - 1/2*t(mu) %*% ginv(s) %*% mu + log(pi_c))
  
}

# calc error using normal computing
lda_calc_error <- function(X, y, vec_mu, vec_pi, S){
  num_class = nlevels(as.factor(y))
  acc = matrix(rep(0, length(y)*num_class), ncol = num_class)
  
  for(i in 1:num_class){
    discr = apply(X, 1, function(x){
      return(lda_linear_discr(x, S, vec_mu[i,], vec_pi[i]))
    })
    acc[,i] = discr
  }
  
  y_hat = apply(acc, 1, which.max)
  
  return(mean((y_hat - min(y_hat)) == (y - min(y))))
  
}

y = convert_class_vector(y)
vec_pi = calc_prior_portions(y)
vec_mu = calc_vec_mu(X, y)
common_s = calc_common_s(X, y, vec_mu)

# calc trainning error
(lda_calc_error(X, y, vec_mu, vec_pi, common_s))

# calc test error
(lda_calc_error(test_data$X, convert_class_vector(test_data$y), vec_mu, vec_pi, common_s))

# using whitening computing
y = convert_class_vector(y)
w_trans = whitenning_transform(X, common_s)
X_star = w_trans$X_w
trans = w_trans$mat_trans
c_s = diag(x=1,nrow = nrow(common_s), ncol = ncol(common_s))
vec_pi = calc_prior_portions(y)
vec_mu = calc_vec_mu(X_star, y)

# Plot the whitening data

plot(X_star[,1], X_star[,2], col=as.numeric(y)+1, xlab = "X1", ylab = "X2", bg=as.numeric(y)+1, pch=21)

# calc trainning error
(lda_calc_error(X_star, y, vec_mu, vec_pi, c_s))

# calc test error
(lda_calc_error(test_data$X %*% trans, convert_class_vector(test_data$y), vec_mu, vec_pi, c_s))


lda_calc_error_wrapper <- function(X_train, y_train, X_test, y_test){
  y = convert_class_vector(y_train)
  vec_mu = calc_vec_mu(X_train, y)
  common_s = calc_common_s(X_train, y, vec_mu)
  
  w_trans = whitenning_transform(X_train, common_s)
  X_star = w_trans$X_w
  trans = w_trans$mat_trans
  c_s = diag(x=1, nrow = nrow(common_s), ncol = ncol(common_s))
  vec_pi = calc_prior_portions(y)
  vec_mu = calc_vec_mu(X_star, y)
  
  # calc training error
  training_error = lda_calc_error(X_star, y, vec_mu, vec_pi, c_s)
  
  # calc test error
  test_error = lda_calc_error(X_test %*% trans, convert_class_vector(y_test), vec_mu, vec_pi, c_s)
  
  return(c(training_error, test_error))
}

# lda with R
df_test = data.frame(test_data$X, convert_class_vector(test_data$y))
names(df_test)[ncol(df_test)] = 'y'
df_train = data.frame(X, y)
lda.fit = lda(y~., data = df_train)
lda.pred = predict(lda.fit, df_train)
mean(lda.pred$class == y)

lda.pred = predict(lda.fit, df_test)
mean(lda.pred$class == df_test$y)
