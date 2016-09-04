library(MASS)

source('gen-data.R')
source('utils.R')

train_data = gen_train_data()
test_data = gen_test_data()
X = train_data$X
y = train_data$y 


qda_linear_discr <- function(x, s, mu, pi_c){
  return(-1/2*t(x-mu) %*% ginv(s) %*% (x-mu) + log(pi_c) - 1/2*log(det(s)))
  
}

qda_calc_error <- function(X, y, vec_mu, vec_pi, list_S){
  num_class = nlevels(as.factor(y))
  acc = matrix(rep(0, length(y)*num_class), ncol = num_class)
  
  for(i in 1:num_class){
    discr = apply(X, 1, function(x){
      return(qda_linear_discr(x, list_S[[i]], vec_mu[i,], vec_pi[i]))
    })
    acc[,i] = discr
  }
  
  y_hat = apply(acc, 1, which.max)
  
  return(mean((y_hat-min(y_hat)) == (y-min(y))))
  
}

qda_calc_error_wrapper <- function(X_train, y_train, X_test, y_test){
  y = convert_class_vector(y_train)
  vec_pi = calc_prior_portions(y)
  vec_mu = calc_vec_mu(X, y)
  separate_s = calc_separate_s(X, y, vec_mu)
  
  trainning_error = qda_calc_error(X, y, vec_mu, vec_pi, separate_s)
  test_error = qda_calc_error(X_test, convert_class_vector(y_test), vec_mu, vec_pi, separate_s)
  
  return(c(training_error, test_error))
}

y = convert_class_vector(y)
vec_pi = calc_prior_portions(y)
vec_mu = calc_vec_mu(X, y)
separate_s = calc_separate_s(X, y, vec_mu)

# calc trainning error
(qda_calc_error(X, y, vec_mu, vec_pi, separate_s))

# calc test error
(qda_calc_error(test_data$X, convert_class_vector(test_data$y), vec_mu, vec_pi, separate_s))

# qda with R
df_test = data.frame(test_data$X, convert_class_vector(test_data$y))
names(df_test)[ncol(df_test)] = 'y'
df_train = data.frame(X, y)
qda.fit = qda(y~., data = df_train)
qda.pred = predict(qda.fit, df_train)
mean(qda.pred$class == y)

qda.pred = predict(qda.fit, df_test)
mean(qda.pred$class == df_test$y)