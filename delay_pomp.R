loc <- url("http://kingaa.github.io/pomp/vignettes/parus.csv")
dat <- read.csv(loc)
head(dat)

library(pomp)
parus <- pomp(dat,times="year",t0=1959)


stochStep <- Csnippet("
  double tx1, tx2, tx3 ;
  tx1 = rnorm(x1,1);
  tx2 = rnorm(tx1,1);
  tx3 = rpois(tx2);

  x1 = tx1; x2 = tx2; x3 = tx3;
  
                      ")

pomp(
  parus,
  rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1),
  initializer=Csnippet("x1 = 1; x2 = 100; x3= 500;"),
  paramnames=c(),
  statenames=c("x1","x2", "x3")
) -> parus

sim <- simulate(parus, params=c(N_0=1,r=12,c=1,sigma=0.5),
                as.data.frame=TRUE, states=TRUE)

plot(x3~time,data=sim,type='o')

rmeas <- Csnippet("pop = rbinom(x3,.5);")

dmeas <- Csnippet("lik = dbinom(pop,x3,q_t,TRUE);")


pomp(parus,
     rmeasure=rmeas,
     dmeasure=dmeas,
     statenames=c("x1","x2","x3"),
     paramnames=c("q_t")
) -> parus


coef(parus) <- c(q_t = .5)


thetastart <- c(q_t=.5)
set.seed(635363)
pmcmc1<-pmcmc(parus, start = thetastart,
              Nmcmc = 5000, Np = 500, max.fail = Inf,
              proposal=mvn.diag.rw(c(q_t = .5)))