normalize_prostate_data <- function(raw_data){
  training_data = raw_data[train == T, ]
  test_data = raw_data[train == F, ]
  
  # extract response and create data frame
  predictors = scale(training_data[,1:8])
  train_center = attr(predictors, 'scaled:center')
  train_scale = attr(predictors, 'scaled:scale')
  train_response_center = mean(training_data[,9])
  
  # trainning data
  lpsa = training_data[,9] - train_response_center
  df_train = data.frame(cbind(predictors, lpsa))
  
  # test data
  predictors = test_data[,1:8]
  for(i in 1:8){
    predictors[,i] = (predictors[,i] - train_center[i])/train_scale[i]
  }
  lpsa = test_data[,9] - train_response_center
  df_test = data.frame(cbind(predictors, lpsa))
  
  return(list(df_train, df_test))
}