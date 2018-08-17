library(surveillance)


get_hohle_estimate <- function(full_reporting_triangle,po_data,D,NSIM){
  #browser()
  
  # create the partially observed triangle (for the other methods i pass in both the fully reported triangle
  # and also the partially observed triangle)
  reporting_triangle <- cbind(full_reporting_triangle[1:(length(full_reporting_triangle)-D),],po_data)
  
  # add 1 because the hohle method requires the triangle to have at least one case per row 
  reporting_triangle <- reporting_triangle +1
  
  count <- 1
  
  # date of hospitalization 
  dHosp <- list()
  
  #date of report
  dReport <- list()
  
  #create a list of dummy date variables to convert to hohle input format which is 
  
  ### --------|--------
  ###  dHosp[1] | dReport[1]
  ###  dHosp[2] | dReport[2]
  ###.......
  
  
  date <- seq(as.Date("2014-1-01"), as.Date("2017-1-1"), by = "2 week")
  
  #iterate through the triangle and take a case from rep[row][col] to map it to dHosp[row] and dReport[row+col-1]
  for (row in 1:nrow(reporting_triangle)){
    for (col in 1:ncol(reporting_triangle)){
       for (ncase in 1:reporting_triangle[row,col]){
         dHosp[[count]] <-date[row]
         dReport[[count]] <- date[row+col-1]
         count <- count + 1
      }
    }
  }
  
  # create the dataframe
  data_df <- do.call(rbind, Map(data.frame, dHosp=dHosp, dReport=dReport))
  
  
  # time NOW is the last reported case
  t.repTriangle <-tail(dReport,1)[[1]]# as.Date("2010-02-13")
  
  #set of nowcast times, in our example it is the last hospital date - D bi weeks
  when <- seq(tail(dHosp,1)[[1]],length.out=D,by="-2 week")
  nc.control <- list(N.tInf.prior=structure("poisgamma",
                                            mean.lambda=100,var.lambda=3000),
                     nSamples=1e2,N.tInf.max =2*max(reporting_triangle))
  
  
  # call the nowcast
  nc.ddcp <- nowcast(now=t.repTriangle,when=when,
                     dEventCol="dHosp",dReportCol="dReport",
                     data=data_df, 
                     method="lawless", D=D,
                     control=nc.control)
  #browser()
  
  
  # for some reason the nowcast method stores the median estimate in the upperbound slot (see surveillance pkg doc)
  tmp <- attr(nc.ddcp,"upperbound")
 
  
  ## it also returns a daily object with a bunch of NAs (because we are on bi-weeks)
  ## so we grab only the filled in estimates from the days corresponding to bi-weeks
  estimate_ <-c()
  for (i in 1:length(tmp)){
    if (!is.na(tmp[i])){
      estimate_ <- c(estimate_ ,tmp[i])
    }
  }
  
  ## right now we just repeat the median
  result_matrix <- matrix(NA,nrow=D,ncol=NSIM)
  
  ##off by1 correction (to undo our transform above )
  for (i in 1:D){
    result_matrix[i,] <- rep(estimate_[i]-1,NSIM)
  }
  return (result_matrix)
}