loc <- url("http://kingaa.github.io/pomp/vignettes/parus.csv")
dat <- read.csv(loc)
head(dat)
dat$pop <- rep(1,nrow(dat))
library(pomp)
parus <- pomp(dat,times="year",t0=1959)


stochStep <- Csnippet("
  double tx1, tx2, tx3 ;
  tx1 = rnorm(x1,100);
  tx2 = rnorm(tx1,100);
  tx3 = rnorm(tx2,100);

  x1 = tx1; x2 = tx2; x3 = tx3;
  
                      ")

pomp(
  parus,
  rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1),
  initializer=Csnippet("x1 = 1; x2 = 1; x3= 1;"),
  paramnames=c(),
  statenames=c("x1","x2", "x3")
) -> parus




rmeas <- Csnippet("pop = rnorm(x3,100);")

dmeas <- Csnippet("lik = dnorm(pop,x3,100*exp(q_t),TRUE);")


pomp(parus,
     rmeasure=rmeas,
     dmeasure=dmeas,
     statenames=c("x1","x2","x3"),
     paramnames=c("q_t")
) -> parus


#coef(parus) <- c(q_t = .5)

pf <- pfilter(parus,Np=1000,params=c(q_t=.5))