library(MHadaptive)
library(MCMCpack)
library(forecast)
library(matrixStats)
library(MASS)
library(boot)
setwd("/Users/gcgibson/reporting_delay/manuscript_code/")

source('utils.R')
source('delay_model.R')
source('zero_model.R')
source('bayes_model.R')


## SET GLOBAL DELAY
D <- 10
## read in data
#reporting_triangle <- read.csv("bangkok_10.csv",header = FALSE)
#reporting_triangle <- generate_data(D)
reporting_triangle <- read.csv("chiang_mai_10.csv",header = FALSE)


### CV START /STOP
start <- 40
stop <-  dim(reporting_triangle)[1]
NSIM <- 100
### INITIALIZE EMPTY VECS
mse_vec <- matrix(NA,nrow=stop-start + 1,ncol=5)
coverage_prob <- matrix(NA,nrow=stop-start + 1,ncol=5)


###  cv cutoffs
for (cv_cutoff in start:stop){
  print (cv_cutoff)
  test_reporting_triangle <- reporting_triangle[1:cv_cutoff,]
  
  ## GET PARTIALLY OBSERVED 
  po_data <- get_po_data(test_reporting_triangle,D,cv_cutoff)
  
  ### MODEL ESTIMATES
  delay_model_estimate <- get_delay_model(test_reporting_triangle,po_data,D,NSIM)
  #zero_inf_offset_2 <- get_zero_model(test_reporting_triangle,po_data,D,cv_cutoff,NSIM,1)
  bayes_us <- get_bayes_model(test_reporting_triangle,po_data,D,cv_cutoff,NSIM,delay_model_estimate)
  bayes_hohle <- get_bayes_model_naive(test_reporting_triangle,po_data,D,cv_cutoff,NSIM,delay_model_estimate)
  
  
  ### COMPUTE MSE
  truth <- rowSums(test_reporting_triangle)[(cv_cutoff-D +1):cv_cutoff]
  
  mse_vec[(cv_cutoff-start +1),1] <- mean((rowMeans(delay_model_estimate)-truth)^2)
  mse_vec[(cv_cutoff-start +1),2] <- mean((rowMeans(bayes_us)-truth)^2)
  
  mse_vec[(cv_cutoff-start +1),3] <- mean((rowMeans(bayes_hohle)-truth)^2)
  #mse_vec[(cv_cutoff-start +1),3] <- mean((rowMeans(zero_inf_offset_2)-truth)^2)
  
 #  png(paste('fcast',toString(cv_cutoff),'.png',sep=""))
 #  plot((cv_cutoff-D+1):cv_cutoff,truth,col="black",type="l",ylim=c(0,max(truth)+1))
 #  lines((cv_cutoff-D+1):cv_cutoff,rowSums(po_data),col="blue")
 #  
 #  #lines((cv_cutoff-D+1):cv_cutoff,rowMeans(zero_inf_offset_2),col="blue")
 #  lines((cv_cutoff-D+1):cv_cutoff,rowMeans(bayes_hohle),col="orange")
 #  
 # lines((cv_cutoff-D+1):cv_cutoff,rowMeans(delay_model_estimate),type="l",col="red")
 #  legend("topleft",legend=c("Truth","Partial", "Bayes","Delay"),
 #         col=c("black","blue", "orange", "red"),lty=c(rep(1,4)), cex=.4)
 #   dev.off()
 #  
  ### Coverage Probability
  quant_delay <-rowQuantiles(delay_model_estimate,probs = c(.025,.975))
  delay_cp <- 0 
  for (cp_iter in 1:D){
    if (quant_delay[cp_iter,1] <= truth[cp_iter] & truth[cp_iter] <=quant_delay[cp_iter,2]){
      delay_cp <- delay_cp + 1
    }
  }
  
  coverage_prob[(cv_cutoff-start +1),1] <- delay_cp/D
  
  
  
  bayes_us_delay <-rowQuantiles(bayes_us,probs = c(.025,.975))
  bayes_us_cp <- 0 
  for (cp_iter in 1:D){
    if (bayes_us_delay[cp_iter,1] <= truth[cp_iter] & truth[cp_iter] <=bayes_us_delay[cp_iter,2]){
      bayes_us_cp <- bayes_us_cp + 1
    }
  }
  coverage_prob[(cv_cutoff-start +1),2] <- bayes_us_cp/D
  
  bayes_hohle_delay <-rowQuantiles(bayes_hohle,probs = c(.025,.975))
  bayes_hohle_cp <- 0 
  for (cp_iter in 1:D){
    if (bayes_hohle_delay[cp_iter,1] <= truth[cp_iter] & truth[cp_iter] <=bayes_hohle_delay[cp_iter,2]){
      bayes_hohle_cp <- bayes_hohle_cp + 1
    }
  }
  coverage_prob[(cv_cutoff-start +1),3] <- bayes_hohle_cp/D
  
  
 
  
  
}

write.matrix(mse_vec,"bangkok_mse")
write.matrix(coverage_prob,"bangkok_coverage_prob")
print(xtable(mse_vec[,1:3], type = "latex"), file = "bangkok.tex")
