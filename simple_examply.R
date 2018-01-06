library(pomp)

loc <- url("https://kingaa.github.io/sbied/intro/parus.csv")
dat <- read.csv(loc)
dat$pop = rep(1,nrow(dat))

parus <- pomp(dat,times="year",t0=1959)


step.fun <- Csnippet("
  double tx1, tx2, tx3;
  tx1 = rnorm(x1,1);
  tx2 = rnorm(tx1,1);
  tx3 = rpois(exp(tx2));

  x1 = tx1; x2 = tx2; x3 = tx3;
")

parus <- pomp(data=dat,time="year",t0=1959,
              rprocess=euler.sim(step.fun=step.fun,delta.t=1/365),
              statenames=c("x1","x2","x3"),paramnames=c("q_1","q_2","q_3"))

rmeas <- Csnippet("
  pop = rbinom(x3,q_1+q_2);")

parus <- pomp(parus,rmeasure=rmeas,statenames="x3",paramnames=c("q_1","q_2","q_3"))

dmeas <- Csnippet("
  lik = dbinom(pop,x3,q_1+q_2,give_log);
                  ")

parus <- pomp(parus,dmeasure=dmeas,statenames="x3",paramnames=c("q_1","q_2","q_3"))


pf <- pfilter(parus,Np=1000,params=c(q_1=.1,q_2=.3,q_3 = .5,x1.0=1,x2.0=1,x3.0=1))

coef(parus) <- c(q_1=.1,q_2=.3,q_3 = .5,x1.0=1,x2.0=1,x3.0=1)


pmcmc(
  pomp(parus,dprior=Csnippet("
                           lik = ddirichlet({q_1,q_2,q_3},1);
                           lik = (give_log) ? lik : exp(lik);"),
       paramnames=c("q_t")),
  Nmcmc=2000,Np=500,verbose=TRUE,
  proposal=mvn.rw.adaptive(rw.sd=c(q_t=1),
                           scale.start=200,shape.start=100)) -> chain



library(coda)
trace <- window(conv.rec(chain,c("q_t")),start=10)
rejectionRate(trace)
effectiveSize(trace)
autocorr.diag(trace)

summary(trace)
plot(trace)

