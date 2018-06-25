

get_bayes_model <- function(test_reporting_triangle,po_data,D,cv_cutoff,NSIM,delay_model_estimate){
  #### mcmc model2 : first get predictions for full nowcasting 
  ### region using arima
  bayes_hohle <- matrix(NA, nrow=D,ncol=NSIM)
  for (s in 1:(NSIM)){
    fcast_training_data <- c(rowSums(test_reporting_triangle)[1:(cv_cutoff-D)],delay_model_estimate[,s][1:(D-offset)])
    fit2 <- auto.arima(fcast_training_data)
    pred_obj <- forecast(fit2,h = offset)$mean
    bayes_estimate_trunc <- c()
    p_star <-   get_p(test_reporting_triangle,D )
    #rdirichlet(1,alpha_star)
    
    
    count <- 2
    for (i in 1:2){
      li_func <- function(n_t_inf,partially_observed, fcast,sum_p){
        tmp <- dbinom(partially_observed,round(n_t_inf),sum_p,log = T) + 
          dnorm(n_t_inf,fcast,100,log=T)
        
        if (!is.finite(tmp)){
          tmp <- -1e100
        }
        return (tmp)
      }
      res_obj <- Metro_Hastings(li_func, c(colSums(po_data)[i+D-2]), prop_sigma = NULL,
                                par_names = NULL, quiet = TRUE,
                                iterations = 1000, burn_in = 100,
                                adapt_par = c(100, 20, 0.5, 0.75),partially_observed=colSums(po_data)[i+D-2],
                                fcast = pred_obj[i],sum_p = sum(p_star[1:count]))
      bayes_estimate_trunc <- c(bayes_estimate_trunc,sample(res_obj$trace,1))
      
      count <- count -1
    }
    
    bayes_hohle[,s] <- c(delay_model_estimate[,s][1:(D-2)],bayes_estimate_trunc)
    
  }
  return (bayes_hohle)
}