#' Bisquare Basis Functions
#'
#' Provides an implementation of the bisquare basis functions defined as
#' \deqn{ \Psi_{j}(s) =
#'  \{ 1 - (||s - c_j||/w)^2\}^2 \mathrm{I}( ||s-c_j|| \leq w ).}
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
#' @return A matrix of bisquare functions evaluated at all combinations of
#'  crd and knots \{ nCrd x r \}.
#'
#' @name bisquare
#' @rdname bisquare
#'
#' @export bisquare
#'
#' @references Cressie, N. and Johannesson, G. (2008). Fixed rank kriging for
#' very large spatial data sets. Journal of the Royal Statistical Society, 
#' Series B, 70, 209--226. <doi:10.1111/j.1467-9868.2007.00633.x>.
#'
#' @examples
#'
#'   data(countyExample)
#'
#'   bisquare(crd = sp::coordinates(county), knots = knots)
#' 
#' @importFrom sp spDists
bisquare <- function(crd, knots, w = NULL, ..., longlat = TRUE) {

  # crd must be a 2 column matrix
  if (!is.matrix(x = crd)) stop("crd must be a matrix")
  if (ncol(x = crd) != 2L) stop("crd must be a matrix with 2 columns")

  # knots must be a 2 column matrix
  if (!is.matrix(x = knots)) stop("knots must be a matrix")
  if (ncol(x = knots) != 2L) stop("knots must be a matrix with 2 columns")

  # if the scaling term is not provided, set as default -- 1.5 * minimum
  # pairwise distance between knots
  if (is.null(x = w)) {
    dis <- sp::spDists(x = knots, y = knots, longlat = longlat)
    w <- 1.5 * min(dis[upper.tri(x = dis, diag = FALSE)])
    message("basis function scale set to ", round(x = w, digits = 4L))
  }

  if (w < 1e-8) stop("w must be positive")

  return( .bisquare(crd = crd, knots = knots, w = w, longlat = longlat) )

}

#' @importFrom sp spDists
.bisquare <- function(..., crd, knots, w, longlat = TRUE) {

  # ||s - cj|| { nCrd x r }
  dis <- sp::spDists(x = crd, y = knots, longlat = longlat)

  return( {1.0 - {dis/w}^2}^2 * {dis <= w} )

}
