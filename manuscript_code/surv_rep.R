data("husO104Hosp")
library(surveillance)
D <- 10
## read in data
reporting_triangle <- read.csv("bangkok_10.csv",header = FALSE)
reporting_triangle <- round(reporting_triangle[1:40,]/1)
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







t.repTriangle <- as.Date("2010-02-13")

when <- seq(t.repTriangle-1,length.out=10,by="-1 day")
nc.control <- list(N.tInf.prior=structure("poisgamma",
                                          mean.lambda=50,var.lambda=3000),
                   nSamples=1e2,N.tInf.max =2*max(reporting_triangle))
nc.ddcp <- nowcast(now=t.repTriangle,when=when,
                   dEventCol="dHosp",dReportCol="dReport",
                   data=data_df, 
                   method="bayes.trunc", D=15,
                   control=nc.control)

#Show time series and posterior median forecast/nowcast
plot(nc,xaxis.tickFreq=list("%d"=atChange,"%m"=atChange),
     xaxis.labelFreq=list("%d"=at2ndChange),xaxis.labelFormat="%d-%b",
     xlab="Time (days)",lty=c(1,1,1,1),lwd=c(1,1,2))



