library(surveillance)
library(MCMCpack)
generate_data <- function(D){
  ts <- c( 173,203,191,145,100,131,81,90,51,101,69,103,175,191,206,166,192,228,259,270,277,286,381,425,381,28,287,279,252,201,209,194,142,207,269,430,566,502,505,676,774,860,1188,1371,1220,1449,1558,2569,4281,5824,4307,133,1169,754,451,375,470,261,119,111,75,54,40,55,67,95,210,247,285,304,291,274,308,358,650,556,4,0)
  
  
  ts_p <- ts[1:26]
  salmonellaDisProg <- create.disProg(week = 1:length(ts_p), 
                                      observed = ts_p,
                                      start = c(1990, 1))
  
  # look at disProg object
  
  
  
  meningo <- disProg2sts(salmonellaDisProg)
  
  
  
  # fit model
  fit <- hhh4(meningo, control = list(ar = list(f = ~ 1),
                                      end = list(f = addSeason2formula(S = 1, period = 26)),
                                      family = "NegBin1"))
  plot(fit)
  
  # simulate from model
  simData <- simulate(fit,nsim=2,simplify = TRUE)
  simDataVec <- c(simData[,,1],simData[,,2])
  data_mat <- matrix(NA,nrow=length(simDataVec),ncol=D)
  p <- rdirichlet(1,c(0.000001,rep(10,5),rep(1,D-6)))
  
  for (row in 1:nrow(data_mat)){
    data_mat[row,] <- p*simData[row]
  }
  #for (i in length(simData):(length(simData)-D)){
    
  #}
  return (data_mat)
  
}



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