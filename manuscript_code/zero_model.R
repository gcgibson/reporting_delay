

get_zero_model <- function(test_reporting_triangle,po_data,D,cv_cutoff,NSIM,offset){
  ###zero inflated model

  ### predicting over last 2 zeros 
  zero_inf <- matrix(NA, nrow=D,ncol=NSIM)
  for (s in 1:NSIM){
    proc_training_data <- c(rowSums(test_reporting_triangle)[1:(cv_cutoff-D)],delay_model_estimate[,s][1:(D-offset)])
    fit <- auto.arima(proc_training_data)
    fcast_l <- forecast(fit,h = offset)
    means <- fcast_l$mean
    upper <- fcast_l$upper
    lower <- fcast_l$lower
    a <- sqrt((upper-means)^2 + (lower -means)^2)
    #a[,2]
    pred_obj <- rnorm(offset,means,a[,2])
    if (offset < D){
      zero_inf[,s] <-c(delay_model_estimate[,s][1:(D-offset)],pred_obj)
    }
    else{
      zero_inf[,s] <-  pred_obj
    }
  }
  return (zero_inf)
}