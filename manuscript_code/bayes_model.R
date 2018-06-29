

get_bayes_model <- function(test_reporting_triangle,po_data,D,cv_cutoff,NSIM,delay_model_estimate_l){
    #### mcmc model2 : first get predictions for full nowcasting 
    ### region using arima
    rep_fact <-200
    
    bayes_hohle <- matrix(NA, nrow=D,ncol=(rep_fact*NSIM))
    

    p_star <-   get_p(test_reporting_triangle,D )
    #rdirichlet(1,alpha_star)
    fcast_training_data <- rowSums(test_reporting_triangle)[1:(cv_cutoff-D)]
    bayes_estimate_trunc <- matrix(NA,nrow=D,ncol=rep_fact)
    
    ## INITIALIZE REVERSE INDEX
    count <- D
    li_func <- function(n_t_inf,partially_observed, fcast,sum_p){
      tmp <- dnorm(partially_observed,round(n_t_inf)*sum_p,10*sqrt(round(n_t_inf)*sum_p*(1-sum_p)),log = T)  + dnorm(n_t_inf,fcast,1000,log=T)
      if (!is.finite(tmp)){
        tmp <- -1e100
      }
      return (tmp)
    }
    
    fit2 <- auto.arima(fcast_training_data)
    pred_obj <- forecast(fit2,h = 1)$mean
    
    res_obj <- Metro_Hastings(li_func, c(rowSums(po_data)[1]), prop_sigma = 5000,
                              par_names = NULL, quiet = TRUE,
                              iterations = 5000, burn_in = 1000,
                              adapt_par = c(500, 20, 0.5, 0.75),partially_observed=rowSums(po_data)[1],
                              fcast = pred_obj[1],sum_p = sum(p_star[1:count])-.000000001)
    bayes_estimate_trunc[1,] <- sample(res_obj$trace,size = rep_fact,replace = TRUE)
    
    count <- count -1
    
    for (i in 2:D){
      ### DEFINE LIKELIHOOD FUNCTION    
      
      
      ##### ITERATED PREDICTIONS
      if (i > 1){
        fcast_training_data <- c(fcast_training_data,sample(bayes_estimate_trunc)[1:i],1)
      }
      #### GET FORECASTS
      fit2 <- auto.arima(fcast_training_data)
      pred_obj <- forecast(fit2,h = 1)$mean
      
      ##### RUN MCMC
      res_obj <- Metro_Hastings(li_func, c(rowSums(po_data)[i]), prop_sigma = 5000,
                                par_names = NULL, quiet = TRUE,
                                iterations = 5000, burn_in = 1000,
                                adapt_par = c(500, 20, 0.5, 0.75),partially_observed=rowSums(po_data)[i],
                                fcast = pred_obj[1],sum_p = sum(p_star[1:count])-.000000001)
      bayes_estimate_trunc[i,] <- sample(res_obj$trace,size = rep_fact,replace = TRUE)
      ##### INVERSE INDEX TO P_SUM
      count <- count -1
    }
    
   
  
  return (bayes_estimate_trunc)
}