

get_bayes_model <- function(test_reporting_triangle,po_data,D,cv_cutoff,NSIM,delay_model_estimate_l,offset){
  #### mcmc model2 : first get predictions for full nowcasting 
  ### region using arima
  rep_fact <-10
  bayes_hohle <- matrix(NA, nrow=D,ncol=(rep_fact*NSIM))
  
  #delay_model_estimate_l<- cbind(delay_model_estimate_l,delay_model_estimate_l)
  for (s in 1:(rep_fact*NSIM)){
    if (s < NSIM){
      delay_index <- s 
    }
    else{
      delay_index <- (s %% NSIM) + 1
    }
    fcast_training_data <- c(rowSums(test_reporting_triangle)[1:(cv_cutoff-D)],delay_model_estimate_l[,delay_index][1:(D-offset)])
    fit2 <- auto.arima(fcast_training_data)
    pred_obj <- forecast(fit2,h = offset)$mean
    bayes_estimate_trunc <- c()
    p_star <-   get_p(test_reporting_triangle,D )
    #rdirichlet(1,alpha_star)
    
    
    count <- offset
    for (i in 1:offset){
      li_func <- function(n_t_inf,partially_observed, fcast,sum_p){
        tmp <- dbinom(partially_observed,round(n_t_inf),sum_p,log = T) + 
          dnorm(n_t_inf,fcast,100,log=T)
        
        if (!is.finite(tmp)){
          tmp <- -1e100
        }
        return (tmp)
      }
      
      
      
      res_obj <- Metro_Hastings(li_func, c(colSums(po_data)[i+D-offset]), prop_sigma = NULL,
                                par_names = NULL, quiet = TRUE,
                                iterations = 1000, burn_in = 100,
                                adapt_par = c(100, 20, 0.5, 0.75),partially_observed=colSums(po_data)[i+D-offset],
                                fcast = pred_obj[i],sum_p = sum(p_star[1:count]))
      bayes_estimate_trunc <- c(bayes_estimate_trunc,sample(res_obj$trace,1))
      
      count <- count -1
    }
    
    if (offset < D){
      bayes_hohle[,s] <- c(delay_model_estimate_l[,delay_index][1:(D-offset)],bayes_estimate_trunc)
    }
    else{
      bayes_hohle[,s] <- bayes_estimate_trunc
    }
  }
  return (bayes_hohle)
}