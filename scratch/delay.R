library(ggplot2)
library(Metrics)
library(MCMCpack)

require(rbiips)

## custom factorization fo joint density
dMN_dim <- function(x,lam,n_t_inf,a1,a2) {
  # Check dimensions of the input and return dimension of the output of
  # distribution dMN
  5
}
dMN_sample <- function(x,lam,n_t_inf,a1,a2,a3) {
  # Draw a sample of distribution dMN
  x1 <- rnorm(1,x,1)
  lam1 <- rnorm(1,x1,1)
  n_t_inf1 <- rpois(1,exp(lam1))
  qt <- rdirichlet(1,c(a1,a2))
  return (c(x1,lam1,n_t_inf1,qt))
  
}
biips_add_distribution('dcustom', 5, dMN_dim, dMN_sample)


# custom dirichlet distribution
ddirichlet_dim <- function(a1,a2,a3) {
  # Check dimensions of the input and return dimension of the output of
  # distribution dMN
  3
}
ddirichlet_sample <- function(a1,a2,a3) {
  # Draw a sample of distribution dMN
  return (rdirichlet(1,c(a1,a2,a3)))
  
}
biips_add_distribution('ddirch', 3, ddirichlet_dim, ddirichlet_sample)




## define data, constants, and initial values  
## here we take T = 4  and D = 3 for illustration so the 
## reporting trapezoid should look like
## n_{4,0}, 0      , 0      , 0
## n_{3,0}, n_{3,1}, 0      , 0
## n_{2,0}, n_{2,1}, n_{2,2}, 0
## n_{1,0}, n_{1,1}, n_{1,2}, 0
T_max = 4
d = 3
n_t_d = matrix(c(
  3, 0,0,0,
  2,2,0,0,
  4,3,4,0,
  1,2,3,0
), nrow =T_max,ncol=T_max,byrow = TRUE)

## compute N_{t,Ts}

N_t_T = c(
  sum(n_t_d[4,1:3]),sum(n_t_d[3,1:3]),sum(n_t_d[2,1:3]),sum(n_t_d[1,1:3])
)






model_file = '/Users/gcgibson/reporting_delay/delay.bug' # BUGS model filename
cat(readLines(model_file), sep = "\n")


t_max = length(N_t_T)
data = list(t_max=t_max, y = N_t_T,mean_x_init=c(1,1,1,1,1))
sample_data = FALSE # Boolean
model = biips_model(model_file, data, sample_data=sample_data) # Create Biips model and sample data
n_part = 1000 # Number of particles
variables = c('x') # Variables to be monitored
mn_type = 'fs'; rs_type = 'stratified'; rs_thres = .5 # Optional parameters
out_smc = biips_smc_samples(model, variables, n_part,
                            type=mn_type, rs_type=rs_type, rs_thres=rs_thres)
diag_smc = biips_diagnosis(out_smc)
summ_smc = biips_summary(out_smc, probs=c(.025, .975))
print (summ_smc$x$f$mean)
