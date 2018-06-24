library(MHadaptive)
library(MCMCpack)
library(forecast)
## read in data
setwd("/Users/gcgibson/reporting_delay/")
reporting_triangle <- read.csv("foo.csv",header = FALSE)
plot(rowSums(reporting_triangle))
D <- 26


start <- 40
stop <- 45
mse_vec <- matrix(NA,nrow=stop-start + 1,ncol=4)
###  cv cutoffs
for (cv_cutoff in start:stop){
  test_reporting_triangle <- reporting_triangle[1:cv_cutoff,]
  #print (test_reporting_triangle)
  
  ##compute dirichlet posterior params
  alpha_star <- colSums(test_reporting_triangle) +1
  
  ## get partially observed data
  po_data <- matrix(NA,nrow=D,ncol=D)
  po_data_index <- 1
  count <- D
  for (i in (cv_cutoff-D+1):cv_cutoff){
    for (j in 1:count){
      po_data[po_data_index,j] <-  test_reporting_triangle[i,j]
    }
    
    for (j in count:D){
      po_data[po_data_index,j] <- 0
    }
    
    
    po_data_index <- po_data_index + 1
    count <- count -1
  }
  
  ## we now have partially observed data
  ## we are ready to start we first get the delay model estimate
  NSIM <- 100
  delay_model_estimate <- matrix(NA,nrow=D,ncol=100)
  count <- D
  for (i in 1:D){
    for (s in 1:NSIM){
      p_star <- rdirichlet(1,alpha_star)
      delay_model_estimate[i,s] <- sum(po_data[i,])/sum(p_star[1:count])
    }
    count <- count -1
  }
  offset <- 2
  ### predicting over last 2 zeros 
  proc_training_data <- rowMeans(delay_model_estimate)[1:(D-offset)]
  fit <- auto.arima(proc_training_data)
  pred_obj <- forecast(fit,h = offset)$mean
  
  ###zero inflated model
  zero_inf <-c(rowMeans(delay_model_estimate)[1:(D-offset)],pred_obj)
  truth <- rowSums(test_reporting_triangle)[(cv_cutoff-D +1):cv_cutoff]
  mse_vec[(cv_cutoff-start +1),1] <- mean((zero_inf-truth)^2)
  mse_vec[(cv_cutoff-start +1),2] <- mean((rowMeans(delay_model_estimate-truth)^2))
  
  #### mcmc model : first get predictions for full nowcasting 
  ### region using arima
  # fcast_training_data <- rowSums(test_reporting_triangle)[1:(cv_cutoff-D)]
  # fit2 <- auto.arima(fcast_training_data)
  # pred_obj <- forecast(fit2,h = D)$mean
  # bayes_estimate <- c()
  # p_star <- rdirichlet(1,alpha_star)
  # 
  # 
  # count <- D
  # for (i in 1:D){
  #   li_func <- function(n_t_inf,partially_observed, fcast,sum_p){
  #     tmp <- dbinom(partially_observed,round(n_t_inf),sum_p,log = T) + 
  #               dnorm(n_t_inf,fcast,1000,log=T)
  #     #print ("---")
  #     #print (n_t_inf)
  #     #print (partially_observed)
  #     #print (fcast)
  #     #print (tmp)
  #     #print (sum_p)
  #     if (!is.finite(tmp)){
  #       tmp <- -1e100
  #     }
  #     return (tmp)
  #   }
  #   res_obj <- Metro_Hastings(li_func, c(colSums(po_data)[i]), prop_sigma = 1000,
  #                  par_names = NULL, iterations = 5000, burn_in = 1000,
  #                  adapt_par = c(100, 20, 0.5, 0.75),partially_observed=colSums(po_data)[i],
  #                  fcast = pred_obj[i],sum_p = sum(p_star[1:count]))
  #   bayes_estimate <- c(bayes_estimate,mean(res_obj$trace))
  # 
  #   count <- count -1
  # }
  
  #### mcmc model2 : first get predictions for full nowcasting 
  ### region using arima
  fcast_training_data <- c(rowSums(test_reporting_triangle)[1:(cv_cutoff-D)],rowMeans(delay_model_estimate)[1:(D-2)])
  fit2 <- auto.arima(fcast_training_data)
  pred_obj <- forecast(fit2,h = 2)$mean
  bayes_estimate_trunc <- c()
  p_star <- rdirichlet(1,alpha_star)
  
  
  count <- 2
  for (i in 1:2){
    li_func <- function(n_t_inf,partially_observed, fcast,sum_p){
      tmp <- dbinom(partially_observed,round(n_t_inf),sum_p,log = T) + 
        dnorm(n_t_inf,fcast,100,log=T)
      print ("---")
      print (n_t_inf)
      print (partially_observed)
      print (fcast)
      print (tmp)
      print (sum_p)
      if (!is.finite(tmp)){
        tmp <- -1e100
      }
      return (tmp)
    }
    res_obj <- Metro_Hastings(li_func, c(colSums(po_data)[i+D-2]), prop_sigma = NULL,
                              par_names = NULL, iterations = 5000, burn_in = 1000,
                              adapt_par = c(100, 20, 0.5, 0.75),partially_observed=colSums(po_data)[i+D-2],
                              fcast = pred_obj[i],sum_p = sum(p_star[1:count]))
    bayes_estimate_trunc <- c(bayes_estimate_trunc,mean(res_obj$trace))
    
    count <- count -1
  }
  
  
  mse_vec[(cv_cutoff-start +1),4] <- mean((c(rowMeans(delay_model_estimate)[1:(D-2)],bayes_estimate_trunc)-truth)^2)
}

