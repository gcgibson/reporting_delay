
### inv logit
inv_transform <- function(x_vec){
  x_vec_exp <- exp(x_vec)
  last_prob = 1/(1+sum(x_vec_exp))
  return(c(x_vec_exp*last_prob,last_prob))
}


get_po_data <- function(test_reporting_triangle,D,cv_cutoff){
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
  return (po_data)
}






#### GET PROBABILITY ESTIMATES FROM TRIANGLE
get_p <- function(test_reporting_triangle,D){
  alpha_star <- colSums(test_reporting_triangle[1:(dim(test_reporting_triangle)[1]-D),]) +1
  
  # transformed_p_hat <- matrix(NA,nrow=dim(test_reporting_triangle)[1],ncol=D)
  # for (i in 1:dim(test_reporting_triangle)[1]){
  #   tmp <- test_reporting_triangle[i,]/sum(test_reporting_triangle[i,])
  #   transformed_p_hat[i,] <- as.numeric(tmp+1)
  # }
  # transformed_p_hat[is.nan(transformed_p_hat)] <- 0
  # 
  # logit_p_s <- matrix(NA,nrow=dim(test_reporting_triangle)[1],ncol=D-1)
  # for (i in 1:dim(test_reporting_triangle)[1]){
  #   tmp <- log(transformed_p_hat[i,1:(D-1)]/transformed_p_hat[i,D])
  #   logit_p_s[i,] <- as.numeric(tmp)
  # }
  # logit_p_s[is.nan(logit_p_s)] <- 0
  #return (inv_transform(mvrnorm(n=1, mu=colMeans(logit_p_s),Sigma = 1*cor(logit_p_s))))
  #new_alpha <- exp(log(alpha_star)+ mvrnorm(n=1,mu=rep(0,D),Sigma = 100*diag(D)))
  return (rdirichlet(1,alpha_star))
}