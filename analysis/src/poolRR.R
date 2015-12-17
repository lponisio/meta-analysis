## takes raw meta-data and a vector of covariates and combines data
## using lajeunesse (2011). Returns a dataframe with keys, the
## observation number, the pooled RR and its variance

poolRR <- function(D, covariates) {
  ## first split the data up into relevant groups:
  if(length(covariates) == 1) {
    dd <- split(D, paste(D$study,
                         D[,covariates],
                         D$year_conv,
                         D$mult_yrs,
                         D$duplicate_trt_dat,
                         sep=';'))
  }
  if(length(covariates) > 1) {
    dd <- split(D, paste(D$study,
                         apply(D[,covariates], 1, paste, collapse=':'),
                         D$year_conv,
                         D$mult_yrs,
                         D$duplicate_trt_dat,
                         sep=';'))
  }

  ## combine response ratios for studies that share treatments using
  ## the method given in Lajeunesse (2011):
  share.trts <- which(sapply(strsplit(names(dd), split=';'),
                             function(x) rev(x)[[1]]) == 'yes')

  dd[share.trts] <- lapply(dd[share.trts],
                           poolCommonControl)

  dd[-share.trts] <-
    lapply(dd[-share.trts], function(x) x[,c('obs', 'lnR', 'varlnR')])
  return(do.call(rbind, dd))
}


