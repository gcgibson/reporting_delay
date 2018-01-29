
library(rstan)

stan_dat <- read_rdump('/Users/gcgibson/reporting_delay/pois_gp.dat.R')
fit_pois <- stan(file="/Users/gcgibson/reporting_delay/pois_gp.stan",   
                 data=stan_dat,
                 iter=200, chains=3)
print(fit_pois, c('rho','alpha','a'))