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


library(ForecastFramework)
library(surveillance)

# load dengue thai data
data_dir <- 'C:/Users/HSDK/Desktop/repos/ForecastFramework/thai-dengue-district-challenge/data/'
dist_counts <- read.csv(paste0(data_dir, 'district-biweek-counts.csv'), header=TRUE)

# get bangkok indices
bangkok_idx <- which(dist_counts$province==10)

# split train/test
n.ahead <- 5
test.idx <- tail(bangkok_idx, n.ahead)
train.data <- dist_counts[bangkok_idx[-test.idx],]
test.data <- dist_counts[bangkok_idx[test.idx],]

############# ############# ############# #############
########## Function to fit HHH4 model ############# #############

fit_hhh4 <- function(train.data, n.ahead, verbose=FALSE){
    # append NA 'hack'
    append.df <- data.frame(province=rep(train.data$province[1], n.ahead),
    district=rep(NA, n.ahead),
    year=rep(tail(train.data$year,1), n.ahead),
    biweek=seq(tail(train.data$biweek,1)+1,
    tail(train.data$biweek,1)+n.ahead, 1) %% 27,
    cases=rep(NA, n.ahead),
    date_sick=rep(NA,n.ahead))
    
    train.data <- rbind(train.data, append.df)
    
    # bangkok case count
    train.disProg <- create.disProg(week=train.data$biweek,
    observed=train.data$cases,
    state=train.data$province,
    start=c(2006,1),
    freq=26,
    neighbourhood=NULL,
    populationFrac=NULL,
    epochAsDate=TRUE)
    
    train.sts <- disProg2sts(train.disProg)
    
    # specify a formula object for the endemic component
    f_S1 <- addSeason2formula(f = ~ 1, S = 1, period = 26)
    
    # fit the Poisson model
    bangkok_NegBinFit <- hhh4(train.sts, control = list(end = list(f = f_S1),
    family = "NegBin1"))
    
    # fit an autoregressive model
    bangkok_NegBinAR <- update(bangkok_NegBinFit, ar = list(f = ~ 1))
    
    # print coefficient estimates
    if(verbose == TRUE) {
        coefs.est <- coef(bangkok_NegBinAR,
        se = TRUE, # also return standard errors
        amplitudeShift = TRUE, # transform sine/cosine coefficients
        idx2Exp = TRUE) # exponentiate remaining parameters
        summary(coefs.est)
    }
    
    # predict ahead
    pred.ahead <- oneStepAhead(bangkok_NegBinAR, nrow(train.data)-n.ahead)
    
}

