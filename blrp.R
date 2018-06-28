library(surveillance)
data(influMen)
# convert to sts class and extract meningococcal disease time series
salmonella <- read.table(system.file("extdata/salmonella.agona.txt", 
                                     package = "surveillance"), header = TRUE)
# look at data.frame
str(salmonella)

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
simData <- simulate(fit,nsim=25)
print (length(simData))

