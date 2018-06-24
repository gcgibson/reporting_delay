
library(MCMCpack)
library(forecast)
## read in data
setwd("/Users/gcgibson/reporting_delay/")
reporting_triangle <- read.csv("foo.csv",header = FALSE)
plot(rowSums(reporting_triangle))
D <- 10


###doMMC

proposalfunction <- function(param){
  return(rnorm(1,mean = param, sd= 10))
}

run_metropolis_MCMC <- function(startvalue, iterations,po_data_l,p_hat,proc_fcast){
  chain = array(dim = c(iterations+1,1))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal,po_data_l,p_hat,proc_fcast) - posterior(chain[i,],po_data_l,p_hat,proc_fcast))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

likelihood <- function(param,po_data_l,p_hat){
  n_t_d_hat = param 
  
  singlelikelihoods = dnorm(po_data_l,p_hat*n_t_d_hat,10000, log = T)
  if (is.nan(singlelikelihoods)){
    singlelikelihoods <- -1e10
    
  }
  return(singlelikelihoods)   
}
prior <- function(param,proc_fcast){
  n_t_d_hat = param 
  n_t_d_hat_prior = dnorm(proc_fcast,sd = 10, log = T)
  return(n_t_d_hat_prior)
}
posterior <- function(param,po_data_l,p_hat,proc_fcast){
  return (likelihood(param,po_data_l,p_hat) + prior(param,proc_fcast))
}
doMCMC <- function(po_data_l,proc_fcast,p_hat){
     
     
      #print (po_data_l)
     
      startvalue = po_data_l
      chain = run_metropolis_MCMC(startvalue, 10000,po_data_l,p_hat,proc_fcast)
      return (chain[100:length(chain)])
}





###  cv cutoffs
for (cv_cutoff in 30:30){
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
  
  proc_training_data <- rowMeans(delay_model_estimate)[1:(D-2)]
  
  fit <- auto.arima(proc_training_data)
  pred_obj <- forecast(fit,h = 2)
  
}

