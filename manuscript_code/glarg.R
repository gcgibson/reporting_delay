library(readr)
#cw<- read_csv("https://raw.githubusercontent.com/FluSightNetwork/cdc-flusight-ensemble/master/model-forecasts/cv-ensemble-models/constant-weights/EW01-2011-constant-weights.csv")
#library(dplyr)
#wk1r1<- filter(cw, location== 'HHS Region 1', target=='1 wk ahead')
values_wk1<- wk1r1$value[1:length(wk1r1$value)-1]

empirical_cdf <- function(q, values){
  bin_index <- 1
  wili <- .1
  while( q > wili){
    bin_index <- bin_index + 1
    wili <- wili +.1
  }
  print(q)
  print (wili)
  print (bin_index)
  return (cumsum(values)[bin_index])
}