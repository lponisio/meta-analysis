source('src/poolRR.R')
source('src/lajeunesse.R')

## makes data arrays for Bayesian analysis; takes the raw
## meta-data set, whether to apply the lajeunesse method to combine RR
## with shared controls, and a vector of explanatory variables. The
## default for covariates is "Crop.species" so RR are never combined
## between crop species. Returns a list of the number of studies
## (Nstudy), the number of years for each study (Nyear), a matrix
## (study, year) of the number of observations in each year (Nobs),
## and an array (study, year, observation) of the RR and their
## variances, and an array of covariates if applicable

makeData <- function(meta.dat,
                     lajeunesse = TRUE,
                     covariates=c('Crop.species')) {

  convert2Int <- function(x) seq_along(x)[match(x, unique(x))]

  if(lajeunesse) {
    dd <- poolRR(meta.dat,
                 covariates=covariates)
    keys.split <- sapply(rownames(dd), strsplit, split=';')
    get <- function(i) as.vector(sapply(keys.split, function(x) x[i]))
    study <- get(1)
    year <- get(3)
  } else {
    dd <- meta.dat
    study <- as.character(meta.dat$study)
    year <- meta.dat$year_conv
  }

  ## make study index
  ind.study <- convert2Int(study)

  ## make year index
  ind.year <- rep(NA, nrow(dd))
  for(i in unique(ind.study))
    ind.year[ind.study==i] <- convert2Int(year[ind.study==i])

  ## make obs index
  ind.obs <- rep(NA, nrow(dd))
  for(i in unique(ind.study))
    for(j in unique(ind.year[ind.study==i]))
      ind.obs[ind.year==j & ind.study==i] <-
        seq_len(sum(ind.year==j & ind.study==i))

  ## make response ratio and vi matrices
  makeMat <- function(d.vec) {
    mat <- array(NA, dim=c(max(ind.study),
                       max(ind.year),
                       max(ind.obs)))
    mat[cbind(ind.study, ind.year, ind.obs)] <- d.vec
    return(mat)
  }
  RR.mat <- makeMat(dd$lnR)
  vi.mat <- makeMat(dd$varlnR)

  ## make index bounds
  nstudy <- max(ind.study)
  nyear  <- apply(RR.mat[,,1],  1, function(x) sum(!is.na(x)))
  nobs   <- apply(RR.mat[,,], 1:2, function(x) sum(!is.na(x)))

  covs <- NA
  num.cats <- NA
  ## covariates (only for lajeunesse)
  if(lajeunesse) {
    num.covs <- length(covariates)
    if(num.covs > 1) {
      covs <- sapply(get(2), strsplit, split=':')
      getCov <- function(cov)
        as.numeric(as.factor(sapply(covs, function(x) x[cov])))
      covs <- lapply(1:num.covs, function(cov)
                     make.mat(getCov(cov)))
      names(covs) <- covariates
      num.cats <- sapply(covs, max, na.rm=TRUE)
    }
  }

  return(list(Nstudy=nstudy,
              Nyear=nyear,
              Nobs=nobs,
              V=vi.mat,
              RR=RR.mat,
              covariates=covs,
              num.cats=num.cats))
}
