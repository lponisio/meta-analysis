
## takes preped data and a scale parameter and packages the data and
## parameters for JAGS, then runs the analysis. Returns the JAGS
## summary output
runAnalysis <- function(dd, scale) {

  my.inits <- function() {
    list()
  }
  my.params <- get.params()

  dd <- list(data=dd, inits=my.inits, params=get.params())

  bugs <- run.model(dd,
                    n.thin=scale,
                    n.iter=(1e4)*scale,
                    n.burnin=1e1*scale,
                    n.chains=3)

  return(list(data=dd,
       bugs=bugs,
       summary=bugs$BUGSoutput$summary))
}
