library(surveillance)


get_hohle_estimate <- function(full_reporting_triangle,po_data,D,NSIM){
  browser()
  reporting_triangle <- cbind(full_reporting_triangle[1:(length(full_reporting_triangle)-D),],po_data)
  reporting_triangle <- reporting_triangle +1
  count <- 1
  dHosp <- list()
  dReport <- list()
  date <- seq(as.Date("2010-1-1"), as.Date("2014-4-8"), by = "day")
  for (row in 1:nrow(reporting_triangle)){
    for (col in 1:ncol(reporting_triangle)){
      if (reporting_triangle[row,col] > 1){
       for (ncase in 1:reporting_triangle[row,col]){
         dHosp[[count]] <-date[row]
         dReport[[count]] <- date[row+col]
         count <- count + 1
       }
      }
    }
  }
  
  data_df <- do.call(rbind, Map(data.frame, dHosp=dHosp, dReport=dReport))
  
  
  
  t.repTriangle <-tail(dReport,1)[[1]]# as.Date("2010-02-13")
  
  when <- seq(tail(dHosp,1)[[1]],length.out=D,by="-1 day")
  nc.control <- list(N.tInf.prior=structure("poisgamma",
                                            mean.lambda=100,var.lambda=3000),
                     nSamples=1e2,N.tInf.max =2*max(reporting_triangle))
  
  
  
  nc.ddcp <- nowcast(now=t.repTriangle,when=when,
                     dEventCol="dHosp",dReportCol="dReport",
                     data=data_df, 
                     method="bayes.trunc", D=D,
                     control=nc.control)
  tmp <- attr(nc.ddcp,"upperbound")#[(length(tmp)-D):length(tmp)]
  # low_ <- tmp[,,1][(length(tmp[,,1])-D):length(tmp[,,1])]
  # high_ <- tmp[,,2][(length(tmp[,,1])-D):length(tmp[,,1])]
  result_matrix <-matrix(NA,nrow=D,ncol=NSIM)
  # estimate_ <-(low_+high_)/2
  for (i in 1:D){
    result_matrix[i,] <- rep(tmp[i],NSIM)
  }
  return (result_matrix)
}