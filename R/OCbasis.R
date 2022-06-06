#' Obled-Cruetin Basis Function
#'
#' Performs a reweighting of radial basis functions ensuring orthonormality.
#'
#' Input 'gbf' allows users to specify a radial basis function beyond
#' the internally implemented bi-square, wendland, and gaussian functions. 
#' If user provides
#' a function, the function must use the following formal arguments: 
#' \itemize{
#' \item crd - the coordinates at which the basis functions are to be evaluated; 
#' \item knots - the knots of the basis functions; 
#' \item w - the scaling factor for the basis function; and
#' \item ... - an ellipsis to avoid argument errors. 
#' }
#' The function must return a matrix of dimension 
#' \{nrow(crd) x nrow(knots)\}.
#' 
#' For completeness, the bi-square functions implemented in the package
#' are of the form
#' \deqn{\Psi_{j}(s) =
#' \{1 - (||s - c_j||/w)^2\}^2 \mathrm{I}( ||s-c_j|| \leq w ).}
#' Note that input 'w' is equivalent to \eqn{w} in the expression above.
#' In addition, if a user were to define an equivalent function 
#' inputs \eqn{s \equiv} 'crd' and \eqn{c_j \equiv} 'knots'.
#'
#' The Wendland basis functions defined as
#' \deqn{ \Psi_{j}(s) =
#'  \{ 1 - d_{j}(s)\}^6 \{35 d_{j}(s)^2 + 18 d_j(s) + 3\}/3 \mathrm{I}( 0 \leq d_{j} \leq 1  ),}
#' where 
#' \deqn{ d_{j}(s) = ||s - c_j||/w.}
#' 
#' The Gaussian radial basis functions defined as
#' \deqn{ \Psi_{j}(s) =
#'  \exp\{- \frac{1}{2} (||s - c_j||/w)^2\}.}
#'
#' @param ... Ignored. Included only to require named inputs.
#'
#' @param spatialData A SpatialPoints object, SpatialPolygons object,
#'   or a list of said objects. The source support data.
#'   If provided as a list, input dB must be integer or SpatialPolygons.
#'
#' @param gbf A function or character object. The function to use to
#'   calculate the radial basis functions of the expansion of the   
#'   Obled-Creutin eigenfunction. The bi-square, wendland, and radial functions
#'   are available through this implementation as 'bisquare', 'wendland' and
#'   'gaussian', respectively. All others must be defined by user.  
#'   See details for further information. 
#'
#' @param dB NULL, integer, or a SpatialPolygons object defining the 
#'   spatial region to be sampled when using Monte Carlo estimates.
#'   If spatialData is a list of spatial objects, dB must be in an integer
#'   specifying the element of spatialData to use as the sampling region
#'   or a SpatialPolygons object.
#'
#' @param w A numeric object. The scaling factor for radial basis functions. 
#'   See details for further information.
#'
#' @param knots A matrix or integer. If a matrix, the knots of the
#'   radial basis functions. If an integer, the number of knots to generate
#'   using fields::cover.design().
#'
#' @param nw An integer object or NULL. The number of MC replicates to 
#'   generate for estimating the O-C eigenfunctions. If <=0 or NULL, 
#'   spatialData must be or include SpatialPoints data.
#'
#' @param nCore An integer object or NULL. The number of cores if parallel 
#'   methods are to be used in the Monte Carlo step.
#'
#' @param longlat A logical object. TRUE if spatialData is 
#'   longitude/latitude data.
#'
#' @return A list containing:
#'  \item{basis}{The \{nSpatial x r\} radial basis function.}
#'  \item{OCnorm}{The \{r x r\} Obled Cruetin weighting matrix.}
#'  \item{knots}{The \{nKnots x 2\} matrix of knots.}
#'  \item{w}{The scaling factor used in basis.}
#'
#' @export
#'
#' @include GBFObj.R obledCruetinBasis.R checkSpatialData.R
#' @include checkNW.R checkDB.R checkInputs.R
#'
#' @name OCbasis
#' @rdname OCbasis
#'
#' @importFrom parallel stopCluster
#'
#' @examples
#' 
#' # create 5x5 square 
#' 
#' poly <- raster::rasterToPolygons(raster::raster(nrows = 5, ncols = 5,
#'                                                 xmn = -1.25, xmx = 1.25,
#'                                                 ymn = -1.25, ymx = 1.25,
#'                                                 res = 0.5,
#'                                                 crs = "+proj=longlat +datum=WGS84"))
#' 
#' df <- data.frame("x" = stats::rnorm(n = 25))
#' 
#' dt <- sp::SpatialPolygonsDataFrame(poly, df)
#' 
#' knots <- expand.grid(c(-0.75,0.0,0.75),c(-0.75,0.0,0.75))
#'
#' OCbasis(spatialData = dt,
#'         gbf = 'bisquare',
#'         knots = knots,
#'         nw = 200L,
#'         nCore = 1L)
#'
#' OCbasis(spatialData = dt,
#'         gbf = 'gaussian',
#'         knots = knots,
#'         nw = 200L,
#'         nCore = 1L)
#'
OCbasis <- function(...,  spatialData,  gbf, knots, dB = NULL, w = NULL,
                    nw = NULL,  nCore = 1L, longlat = TRUE) {

  # ensure that spatialData is a SpatialPointsDataFrame, 
  # SpatialPolygonsDataFrame or a list of said objects
  spatialData <- .checkSpatialData(spatialData = spatialData)

  # ensure appropriate nw value
  nw <- .checkNW(spatialData = spatialData, nw = nw)

  # ensures that dB is one of NULL, integer or SpatialPolygons
  dB <- suppressMessages(expr = .checkDB(dB = dB, 
                                         spatialData = spatialData, 
                                         cage = TRUE))

  # if dB SpatialPolygons, return SpatialPolygons 
  # if dB integer and spatialData is a list, returns integer
  # if dB NULL and spatialData is SpatialPoints or SpatialPolygons, NULL
  # if dB NULL and spatialData is a list, error
  dB <- dB$dB

  # object of class GBFObj containing the weight, knots, and
  # function name for calculating the basis functions.
  gbfObj <- .verifyGBF(gbf = gbf, 
                       weight = w, 
                       knots = knots, 
                       dB = dB,
                       longlat = longlat,
                       spatialData = spatialData)

  # establish cluster if using parallel methods
  localCluster <- .check_nCore(nCore = nCore, parallelLog = NULL)

  # calculate basis functions
  phi <- .obledCruetinBasis(spatialData = spatialData,
                            dB = dB,
                            gbfObj = gbfObj,
                            nw = nw,
                            localCluster = localCluster,
                            verify = TRUE)


  # end cluster if created
  if (!is.null(x = localCluster)) {
    parallel::stopCluster(cl = localCluster)
  }

  return( list("basis" = phi$basis, 
               "OCnorm" = phi$OCnorm,  
               "knots" = phi$gbfObj@args$knots,
               "w" = phi$gbfObj@args$w) )
}
