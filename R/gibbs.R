#  Gibbs Sampling
#
#  Function is not exported
#
#  @param Z A numeric object. The process of interest. { nSource }
#
#  @param X A matrix object. The covariates for the large-scale model. 
#    { nSource x p }
#
#  @param H A matrix object. Indicator of overlap of spatial and dB objects.
#    { nSource x nB }
#
#  @param psi A matrix object. The OC basis. { nSource x r }
#
#  @param sigma2_eps A numeric vector or matrix object. The survey variance. 
#    { nSource } {nSource x nSource}
#
#  @param sigma2_beta A numeric vector object. The sigma^2 of beta 
#    distribution. {p}
#
#  @param sigma2_xi A numeric object. The sigma^2 of xi distribution. {1}
#
#  @param qObj A QPriorObj object. The prior distribution for Q.
#
#  @param gibbsKeep An integer vector object. The Gibbs samples to keep.
#    {nKept}
#
#  @param ffdir A character object or NULL. The directory in which 
#    Gibbs results are stored.
#
# @return a list containing
#  \item beta { p x nGibbs } an (ff) matrix of estimated beta parameters
#  \item eta { r x nGibbs } an (ff) matrix of estimated eta parameters
#  \item xi { n x nGibbs } an (ff) matrix of estimated xi parameters
#
#' @importFrom ff ff
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm
#' @importFrom pracma pinv
#' @import Matrix
#'
#' @include QPriorObj.R
.gibbs <- function(...,
                   Z,
                   X,
                   H,
                   psi,
                   sigma2_eps,
                   sigma2_beta, 
                   sigma2_xi,
                   qObj,
                   gibbsKeep,  
                   ffdir) {

  message("generating Gibbs samples")

  isBad <- function(x, nm) {
             if (any(is.na(x = x)) || 
                 tryCatch(any(is.nan(x = x)),
                          error = function(e){FALSE}) || 
                 any(is.infinite(x = x))) {
               stop(nm, 'contains invalid values (Inf, NA, or NaN)', 
                    call. = FALSE)
             }
           }

  #### confirm appropriate inputs
  isBad(x = Z, nm = 'variable of interest')

  isBad(x = sigma2_eps, nm = 'sigma2_eps')
  isBad(x = 1.0 / sigma2_eps, nm = 'sigma2_eps')

  isBad(x = sigma2_beta, nm = 'sigma2_beta')
  isBad(x = 1.0 / sigma2_beta, nm = 'sigma2_beta')

  isBad(x = sigma2_xi, nm = 'sigma2_xi')
  isBad(x = 1.0 / sigma2_xi, nm = 'sigma2_xi')

  #### dimensions of the problem

  # number of covariates in X {p}
  p <- ncol(x = X)

  # number of radial basis functions in expansion of phiOC {r}
  r <- ncol(x = psi)

  # number of spatial units {nSource}
  n <- length(x = Z)

  # number of Gibbs samples that are kept {nKept}
  nKept <- length(x = gibbsKeep)

  # number of units in minimum resolution {nB}
  nB <- ncol(x = H)

  #### storage matrix for Gibbs
  xi <- .makeStorage(ffdir = ffdir, nrow = nB, ncol = nKept)
  beta <- .makeStorage(ffdir = ffdir, nrow = p, ncol = nKept)
  eta <- .makeStorage(ffdir = ffdir, nrow = r, ncol = nKept)

  #### one time calculations to speed up function

  # {nSource x nSource} matrix object
  if (is.matrix(x = sigma2_eps)) {
    sig2epsInv <- pracma::pinv(A = sigma2_eps)
  } else {
    sig2epsInv <- diag(x = 1.0/sigma2_eps, nrow = n, ncol = n)
  }

  ### beta distribution

  # initialize beta

  tempBeta <- drop(x = tcrossprod(x = pracma::pinv(A = crossprod(x = X)), 
                                  y = X) %*% Z)

  XBeta <- X %*% tempBeta

  tB <- tempBeta / sigma2_beta

  # {p x nSource}
  XpinvV <- t(x = X) %*% sig2epsInv

  # { p x p }
  temp <- diag(x = {1.0 / sigma2_beta}, nrow = p, ncol = p) + XpinvV %*% X
  denoBeta <- tryCatch(expr = pracma::pinv(A = temp),
                       error = function(e) {
                                 stop("unable to obtain beta denominator\n", 
                                      e$message, call. = FALSE)
                               })

  isBad(x = denoBeta, nm = 'denoBeta')

  denoBetaX <- denoBeta %*% XpinvV

  tB <- denoBeta %*% tB

  ### xi distribution

  # initialize xi
  # {nB}
  tempXi <- numeric(length = nB)

  # expand tempXi to source data
  # all points in Bi have the same xi (1s pull the correct one)
  # all polygons containing xi have a weighted value of xi
  #   intersecting area / total area
  # {nSource}
  HXi <- H %*% tempXi

  # {nB x nB} Matrix object
  Zsig2epsInvZ <- Matrix::t(x = H) %*% sig2epsInv %*% H

  # {nB x nB} matrix object
  xiInv <- diag(x = 1.0/sigma2_xi, nrow = nB, ncol = nB)

  # {nB x nB} Matrix object
  U <- xiInv + Zsig2epsInvZ

  # eigen object
  eigenU <- eigen(x = U, symmetric = TRUE)

  # {nB x nSource} Matrix object
  b <- eigenU$vectors %*% diag(x = 1.0/eigenU$values) %*% 
       t(x = eigenU$vectors) %*% Matrix::t(x = H) %*% sig2epsInv

  # {nB x nB} matrix object
  denoXi_sr <- eigenU$vectors %*% diag(x = 1.0/sqrt(x = eigenU$values)) %*% 
               t(x = eigenU$vectors)


  # {nB} Matrix object
  sigmaZ <- b %*% Z
  # {nB x nBeta} Matrix object
  sigmaX <- b %*% X
  # {nB x r} Matrix object
  sigmaS <- b %*% psi

  ### eta distribution

  # initialize eta
  tempEta <- stats::rnorm(n = r, mean = 0.0, sd = 1.0)

  # {nSource x 1} matrix object
  SEta <- psi %*% tempEta

  # { r x nSource } Matrix object
  SpinvV <- Matrix::Matrix(t(x = psi) %*% sig2epsInv)

  # { r x r } Matrix object
  SpinvVS <- SpinvV %*% psi


  #### Sampling
  ind <- 1L

  nGibbs <- gibbsKeep[nKept]

  for (ng in 2L:nGibbs) {

    if (ng %% 500L == 0L) {
      message(ng, " ", appendLF = FALSE)
      if (ng %% 5000L == 0L) message("")
    }

    ### full conditional beta
    # { p }
    mubeta <- denoBetaX %*% {Z - HXi - SEta} + tB

    tempBeta <- MASS::mvrnorm(n = 1L, mu = mubeta, Sigma = denoBeta)

    isBad(x = tempBeta, nm = 'beta')

    XBeta <- X %*% tempBeta

    ### full conditional xi
    # { nB }

    # {nB}
    muxi <- sigmaZ - sigmaX %*% tempBeta - sigmaS %*% tempEta

    # {nB}
    tempXi <- muxi + denoXi_sr %*% stats::rnorm(n = nB, mean = 0.0, sd = 1.0)

    isBad(x = tempXi, nm = 'xi')

    HXi <- H %*% tempXi

    ### full conditional eta
    # { r }

    qObj <- .gibbsQ(qObj = qObj, 
                    SpinvVS = SpinvVS, 
                    SpinvV = SpinvV)

    muEta <- qObj@mu %*% {Z - XBeta - HXi}

    tempEta <- muEta + qObj@deno %*% stats::rnorm(n = r, mean = 0.0, sd = 1.0)

    isBad(x = tempEta, nm = 'eta')

    SEta <- psi %*% tempEta

    qObj <- .metroQ(qObj = qObj, eta = tempEta)

    if (ng == gibbsKeep[ind]) {
      beta[,ind] <- as.vector(x = tempBeta)
      eta[,ind] <- as.vector(x = tempEta)
      xi[,ind] <- as.vector(x = tempXi)
      ind <- ind + 1L
    }

  }

  message("")

  return( list("beta" = beta, "eta" = eta, "xi" = xi) )

}

.makeStorage <- function(ffdir, nrow, ncol) {

  if (is.null(x = ffdir)) {
    return( tryCatch(expr = matrix(data = 0.0, nrow = nrow, ncol = ncol),
                     error = function(e){
                               stop("unable to create matrix of dim", 
                                    nrow, ncol, e$message, call. = FALSE)
                             }) )
  } else {
    return( ff(vmode = "double", 
               dim = c(nrow, ncol),
               dimnames = list(1L:nrow, 1L:ncol)) )
  }
}
