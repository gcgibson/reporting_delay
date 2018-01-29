library(pomp)

pompExample(ou2)

data(ou2)
ou2.dat <- as.data.frame(ou2)
pomp(
  data=ou2.dat[c("time","y1","y2")],
  times="time",
  t0=0,
  rprocess=discrete.time.sim(
    step.fun=function (x, t, params, ...) {
      eps <- rnorm(n=2,mean=0,sd=1) # noise terms
      xnew <- c(
        x1=params["alpha.1"]*x["x1"]+params["alpha.3"]*x["x2"]+
          params["sigma.1"]*eps[1],
        x2=params["alpha.2"]*x["x1"]+params["alpha.4"]*x["x2"]+
          params["sigma.2"]*eps[1]+params["sigma.3"]*eps[2]
      )
      names(xnew) <- c("x1","x2")
      xnew
    }
  )
) -> ou2.Rplug

pmcmc(
  pomp(ou2,dprior=Csnippet("
                           lik = dnorm(q_t,-0.5,1,1);
                           lik = (give_log) ? lik : exp(lik);"),
       paramnames=c("q_t")),
  Nmcmc=20,Np=500,verbose=TRUE,
  proposal=mvn.rw.adaptive(rw.sd=c(q_t=1),
                           scale.start=200,shape.start=100)) -> chain
continue(chain,Nmcmc=20,proposal=mvn.rw(covmat(chain))) -> chain
plot(chain)
chain <- pmcmc(chain)
plot(chain)

library(coda)
trace <- window(conv.rec(chain,c("alpha.2","alpha.3")),start=20)
rejectionRate(trace)
effectiveSize(trace)
autocorr.diag(trace)

summary(trace)
plot(trace)