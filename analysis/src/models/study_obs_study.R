run.model <- function(d,
                      n.burnin,
                      n.iter,
                      n.thin,
                      n.chains) {
  
  sink('model.jags')
  cat('model { 
    for(study in 1:Nstudy) {
      for(year in 1:Nyear[study]) {
        for(obs in 1:Nobs[study,year]) {
          P[study,year,obs] <- 1/V[study,year,obs]
          RR[study,year,obs] ~ dnorm(mu.RR[study,year,obs],
                                     P[study,year,obs])
          mu.RR[study,year,obs] ~ dnorm(mu.study[study],
                                    tau.yr.obs[study])
        }
      }
      mu.study[study] ~ dnorm(mu, tau)
      tau.yr.obs[study] ~ dgamma(shape.obs, scale.obs)
    }

    mu ~ dnorm(0, 1e-4)
    exp.mu <- exp(mu)
    tau <- 1 / (sigma * sigma)
    sigma ~ dunif(0, 100)

    shape.obs <- (1/cv.obs)^2
    cv.obs ~ dunif(0, 100)
    scale.obs <- (1/in.scale.obs)^2
    in.scale.obs ~ dunif(0, 100)
  }', fill = TRUE)
  sink()
  
  attach(d$data)
  res <- jags.parallel(data=list('Nstudy', 'Nyear', 'Nobs', 'V', 'RR'),
                       inits=d$inits,
                       parameters.to.save=d$params,
                       model.file='model.jags',
                       n.burnin=n.burnin,
                       n.chains=n.chains,
                       n.iter=n.iter,
                       n.thin=n.thin,
                       working.directory=NULL)
  detach(d$data)
  
  res
}

## specify the parameters to be monitored
get.params <- function()
  c('exp.mu', 'sigma',
    'cv.obs')
