#' Wendland Basis Functions
#'
#' Provides an implementation of the Wendland basis functions defined as
#' \deqn{ \Psi_{j}(s) =
#'  \{ 1 - d_{j}(s)\}^6 \{35 d_{j}(s)^2 + 18 d_j(s) + 3\}/3 \mathrm{I}( 0 \leq d_{j} \leq 1  ),}
#' where 
#' \deqn{ d_{j}(s) = ||s - c_j||/w.}
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
#' @return A matrix of Wendland functions evaluated at all combinations of
#'  crd and knots \{ nCrd x r \}.
#'
#' @name wendland
#' @rdname wendland
#'
#' @import sp
#'
#' @export wendland
#'
#' @references Wendland, H. (1998). Error estimates for interpolation by
#'   compactly supported radial basis functions of minimal degree. Journal
#'   of Approximation Theory, 93,258-272. <doi:10.1006/jath.1997.3137>.
#'
#' @examples
#'
#'   data(countyExample)
#'
#'   wendland(crd = sp::coordinates(county), knots = knots)
#' 
wendland <- function(crd, knots, w = NULL, ..., longlat = TRUE) {

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

  return( .wendland(crd = crd, knots = knots, w = w, longlat = longlat) )

}

#' @importFrom fields rdist.earth
.wendland <- function(..., crd, knots, w, longlat = TRUE) {

  # ||s - cj|| { nCrd x r }
  if (longlat) {
    dis <- fields::rdist.earth(x1 = crd, x2 = knots, miles = FALSE) / w
  } else {
    dis <- sp::spDists(x = crd, y = knots, longlat = longlat) / w
  }

  nzero <- dis <= 1.0
  dis2 <- dis[nzero]
  tt <- {1.0 - dis2}^6 * {35.0*dis2^2 + 18.0*dis2 + 3.0} / 3.0
  dis[] <- 0.0
  dis[nzero] <- tt

  return( dis )

}
