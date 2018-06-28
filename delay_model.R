


get_delay_model <- function(test_reporting_triangle,po_data,D,NSIM){
  delay_model_estimate <- matrix(NA,nrow=D,ncol=NSIM)
count <- D
  for (i in 1:D){
    for (s in 1:NSIM){
      p_star <- get_p(test_reporting_triangle,D )
      #rdirichlet(1,alpha_star)
      tmp_p <- sum(p_star[1:count])
      
      # logit strategy for adding noise to proportion observed 
      # if (tmp_p < 1 & tmp_p>0){
      #  
      #   tmp_p <- log(tmp_p/(1-tmp_p)) + rnorm(1,0,10)
      #   tmp_p <- inv.logit(tmp_p)
      # }
      #print (tmp_p)
      
      ## normal noise strategy 
      delay_model_estimate[i,s] <- sum(po_data[i,])/tmp_p#rnorm(1,)#/(tmp_p),50)
    }
    count <- count -1
  }
   return (delay_model_estimate)
}