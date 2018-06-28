

get_bayes_model_descending_var <- function(test_reporting_triangle,po_data,D,cv_cutoff,NSIM,delay_model_estimate_l){
  #### mcmc model2 : first get predictions for full nowcasting 
  ### region using arima
  rep_fact <-1000
  offset <- D
  bayes_hohle <- matrix(NA, nrow=D,ncol=(rep_fact*NSIM))
  
  #delay_model_estimate_l<- cbind(delay_model_estimate_l,delay_model_estimate_l)
  for (s in 1:1){
    
    
    p_star <-   get_p(test_reporting_triangle,D )
    #rdirichlet(1,alpha_star)
    fcast_training_data <- rowSums(test_reporting_triangle)[1:(cv_cutoff-D)]
    bayes_estimate_trunc <- matrix(NA,nrow=D,ncol=rep_fact)
    
    count <- offset
    
    
    
    
    for (i in 1:offset){
      li_func <- function(n_t_inf,partially_observed, fcast,sum_p){
        tmp <- dbinom(partially_observed,round(n_t_inf),sum_p,log = T) 
        if (i < D-2){
         tmp <- tmp+ dnorm(n_t_inf,fcast,10000,log=T)
        }
        else{
          tmp <- tmp + dnorm(n_t_inf,fcast,100,log=T)
        }
        
        if (!is.finite(tmp)){
          tmp <- -1e100
        }
        return (tmp)
      }
      
      
      fcast_training_data <- c(fcast_training_data,rowMeans(bayes_estimate_trunc))
      
      fit2 <- auto.arima(fcast_training_data)
      pred_obj <- forecast(fit2,h = offset)$mean
      
      
      res_obj <- Metro_Hastings(li_func, c(rowSums(po_data)[i+D-offset]), prop_sigma = 5000,
                                par_names = NULL, quiet = TRUE,
                                iterations = 5000, burn_in = 1000,
                                adapt_par = c(500, 20, 0.5, 0.75),partially_observed=rowSums(po_data)[i+D-offset],
                                fcast = pred_obj[i],sum_p = sum(p_star[1:count])-.000000001)
      bayes_estimate_trunc[i,] <- sample(res_obj$trace,size = rep_fact,replace = TRUE)
      
      count <- count -1
    }
    
   
  }
  return (bayes_estimate_trunc)
}