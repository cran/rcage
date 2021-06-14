#' Gaussian Radial Basis Functions
#'
#' Provides an implementation of the Gaussian radial basis functions defined as
#' \deqn{ \Psi_{j}(s) =
#'  \exp\{- \frac{1}{2} (||s - c_j||/w)^2\}.}
#'
#' Distances between reference coordinates and knots are obtained using 
#'   sp::spDists().
#'
#' @param crd A matrix object. The (x,y) coordinates of the reference
#'   points \{nCrd x 2\}.
#'
#' @param knots A matrix object. The (x,y) coordinates of the knots \{r x 2\}.
#'
#' @param w A numeric object. The positive scaling factor (bandwidth).
#'
#' @param ... ignored. Included only to require naming of inputs that follow.
#'
#' @param longlat A logical object. If FALSE, Euclidean
#'    distance is calculated; if TRUE, Great Circle distance is calculated.
#'    See ?sp::spDists for more information.
#'
#' @return A matrix of Gaussian functions evaluated at all combinations of
#'  crd and knots \{ nCrd x r \}.
#'
#' @name radial
#' @rdname radial
#'
#' @import sp
#'
#' @export radial
#'
#' @examples
#'
#'   data(countyExample)
#'
#'   radial(crd = sp::coordinates(county), knots = knots)
#' 
radial <- function(crd, knots, w = NULL, ..., longlat = TRUE) {

  if (!is.matrix(x = crd)) stop("crd must be a matrix")
  if (ncol(x = crd) != 2L) stop("crd must be a matrix with 2 columns")

  if (!is.matrix(x = knots)) stop("knots must be a matrix")
  if (ncol(x = knots) != 2L) stop("knots must be a matrix with 2 columns")

  if (is.null(x = w)) {
    dis <- sp::spDists(x = knots, y = knots, longlat = longlat)
    w <- 1.5 * min(dis[upper.tri(x = dis, diag = FALSE)])
    message("basis function scale set to ", round(w,4))
  }

  if (w < 1e-8) stop("w must be positive")

  return( .gauss(crd = crd, knots = knots, w = w, longlat = longlat) )

}

.gauss <- function(..., crd, knots, w, longlat = TRUE) {

  # ||s - cj|| { nCrd x r }
  dis <- sp::spDists(x = crd, y = knots, longlat = longlat) / w

  return( exp(x = -dis^2 * 0.5) )

}
