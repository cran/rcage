# Inverse of the MI Prior Distribution
#
# Calculates the inverse of the MI prior as defined in the supplement to 
#   original manuscript and returns its inverse. See reference for details.
#
# Method is not exported.
#
# For clarity -- function returns 
#   \deqn{(\sigma^2 Q)^{-1} = R_{B}^{-1} A^{+}{Q_B'(I-A)Q_{B}} R_{B}^{-1}}
#
# @references Bradley, J. R., Wikle, C. K., and Holan, S. H. (2017).
#   Regionalization of Multiscale Spatial Processes using a Criterion
#   for Spatial Aggregation Error. Journal of the Royal Statistical Society -
#   Series B, 79, 815--832. <doi:10.1111/rssb.12179>
#
# @param ... Ignored.
#
# @param psi A matrix object. The OC basis. {nB x r}
#
# @param dB A SpatialPolygons or SpatialPoints object. The finest
#   resolution. {nB}
#
# @param dMax A numeric object. Maximum distance between two points to
#   be considered neighbors. Ignored if 'dB' is SpatialPolygons.
#
# @return matrix { r x r }
#
#' @importFrom rgeos gTouches
#' @importFrom pracma pinv
#' @importFrom stats quantile
#
#' @include semiDefinite.R

.miPrior <- function(..., psi, dB, dMax = NULL) {

  # number of areal units in finest resolution
  nB <- nrow(x = psi)

  # number of basis functions in expansion
  r <- ncol(x = psi)

  if (length(x = dB) != nB) {
    stop("contact developer with code miPrior1", call. = FALSE)
  }

  message("calculating MI prior")

  # QR decomposition of phiOC
  qrResult <- tryCatch(expr = qr(x = psi), 
                       error = function(e) { 
                                 stop("unable to obtain QR decomp\n",
                                      e$message, call. = FALSE) 
                               })

  Q <- qr.Q(qr = qrResult)
  R <- qr.R(qr = qrResult)

  addI <- r - ncol(x = Q)

  if (addI > 0L) {

    Q <- cbind(Q, matrix(data = 0.0, nrow = nrow(x = Q), ncol = addI))
    R <- rbind(R, matrix(data = 0.0, nrow = addI, ncol = ncol(x = R)))

  }

  # adjacency matrix for areal units
  if (is(object = dB, class2 = "SpatialPolygons")) {
    # if dB is spatialPolygons object, adjacency is 0/1 matrix indicating
    # shared boundaries
    adjacency <- tryCatch(expr = rgeos::gTouches(spgeom1 = dB, byid = TRUE),
                          error = function(e){
                                    stop("unable to obtain adjacency matrix\n",
                                         e$message, call. = FALSE)
                                  })
  } else {
    # if dB is spatialPoints object, adjacency is 0/1 matrix indicating
    # points within distance of dMax
    D <- sp::spDists(x = dB)
    if (is.null(x = dMax)) {
      dMax <- stats::quantile(x = D, probs = 0.1)
      message("adjacent points are within ", round(dMax, 4))
    }
    adjacency <- D < dMax
  }

  # objects cannot be adjacent to themselves
  diag(x = adjacency) <- FALSE

  adjacency <- unname(obj = adjacency)

  # { r x r }
  qInv <- t(x = Q) %*% {diag(nrow = nB) - adjacency} %*% Q

  # approximate best Hermitian positive semi-definite
  aPlus <- .semiDefinite(A = qInv)

  # inverse of the R-matrix of QR decomposition { r x r }
  ginvR <- tryCatch(expr = pracma::pinv(A = R), 
                    error = function(e) {
                              stop("unable to invert R\n", 
                                   e$message, call. = FALSE)
                            })

  return( ginvR %*% aPlus %*% ginvR )

}
