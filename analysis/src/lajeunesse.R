## applies the Lajeunesse (2011) method for combining RR with shared
## controls. Returns a pooled RR and variance
poolCommonControl <- function(dd) {

  ## function to make covariance matrix:
  covMat <- function(x){
    ## case 1: organic is repeated
    if(all(diff(as.numeric(x[,"mean.org"]))==0))
      vals <- as.numeric(x[1, c("sd.org", "n.organic", "mean.org")])
    if(all(diff(as.numeric(x[,"mean.conv"]))==0))
      vals <- as.numeric(x[1, c("sd.conv", "n.conv", "mean.conv")])

    m <- matrix(vals[1]^2/(vals[2]*vals[3]^2),
                ncol=nrow(dd), nrow=nrow(dd))
    diag(m) <- as.numeric(x[,"varlnR"])
    m
  }

  ## covariance matrix and its inverse:
  V <- covMat(dd)
  V.inv <- solve(V)

  X <- matrix(1, nrow(V))
  E <- matrix(dd[,"lnR"], nrow(V))
  ## compute and return RR_bar and variance(RR_bar), as defined in
  ## Eq. 3 and subsequent text in Lajeunesse (2011):
  return(c(obs=dd$obs[1],
    lnR=as.numeric(solve(t(X) %*% V.inv %*% X) %*%
      (t(X) %*% V.inv %*% E)),
    varlnR=as.numeric(solve(t(X) %*% V.inv %*% X))))
}
