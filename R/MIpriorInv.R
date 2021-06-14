#' Inverse of the MI Prior Distribution
#'
#' Calculates the inverse of the MI prior distribution, which is motivated
#'   by specifying the covariance so that it is "close" to the covariance from 
#'   an Intrinsic Conditional Auto-Regressive model on the domain of the finest
#'   resolution. See reference for details.
#'
#' For clarity -- function returns 
#'   \deqn{K^{-1} = R_{B}^{-1} A^{+}{Q_B'(I-A)Q_{B}} R_{B}^{-1},}
#'  where \eqn{A^{+}} is the first order adjacency matrix and \eqn{Q_B, R_B} is
#'  the QR decomposition of the basis matrix.
#'
#' @references Bradley, J. R., Wikle, C. K., and Holan, S. H. (2017).
#'   Regionalization of Multiscale Spatial Processes using a Criterion
#'   for Spatial Aggregation Error. Journal of the Royal Statistical Society -
#'   Series B, 79, 815--832. <doi:10.1111/rssb.12179>
#'
#' @param psi A matrix object. The estimated OC basis \{nB x r\}, where nB
#'    is the number of areal units in the finest resolution spatial object (dB)
#'
#' @param dB A SpatialPolygons or SpatialPoints object. The finest 
#'    resolution \{nB\}
#'
#' @param dMax Numeric maximum distance between points to be considered
#'   adjacent. Ignored if dB is SpatialPolygons. If dB is SpatialPoints and
#'   dMax is not specified, it is taken to be the 0.1 quantile of the distances.
#'
#' @param ... Ignored.
#'
#' @return The inverse of the MI prior as an rxr matrix.
#'
#' @name MIpriorInv
#' @rdname MIpriorInv
#'
#' @export
#'
#' @include miPrior.R
#'
#' @examples
#'
#' data(countyExample)
#'
#' nc <- county[county@data[,"STATE"] == 37, ]
#' psi <- matrix(data = rbinom(n = 1000, size = 1, prob = 0.5), 
#'               nrow = 100L, ncol = 10L)
#' MIpriorInv(psi = psi, dB = nc)
#'
MIpriorInv <- function(psi, dB, ..., dMax = NULL) {

  if (!is.matrix(x = psi)) stop("psi must be a matrix")

  if (!is(object = dB, class2 = "SpatialPoints") &&
      !is(object = dB, class2 = "SpatialPolygons")) {
    stop("dB must be a SpatialPoints or SpatialPolygons object")
  }

  if (nrow(x = psi) != length(x = dB)) {
    stop("dim error - rows of psi must match length of dB")
  }

  return( .miPrior(psi = psi, dB = dB, dMax = dMax) )

}
