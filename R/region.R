#' Regionalization of Multiscale Spatial Processes
#' 
#' Regionalization of multiscale spatial processes based on a criterion 
#' for spatial aggregation error. The multiscale representation is 
#' a truncated Karhunen-Loaeve expansion using Obled-Creutin 
#' eigenfunctions. The method is incorporated within a Bayesian
#' framework using a Markov chain Monte Carlo implementation of a 
#' latent spatial model. This function differs from 
#' optRegion() in that clustering algorithms are not used
#' to determine an optimal clustering, but spatial regions are clustered
#' according to the provided regionalization and CAGE/DCAGE is calculated.
#' 
#' Input 'spatialData' must be a single spatial object or a list of 
#' spatial objects as defined by the \code{sp}
#' package. If it is a SpatialPointsDataFrame object and nw = 0/NULL, 
#' the package will use the point data to obtain Psi and W. 
#' If nw>0, the package will use Monte Carlo to estimate Psi and W. 
#' 
#' If input 'spatialData' does not represent the finest resolution
#' upon which inference is made, input 'dB' must be set to specify
#' the finest resolution. The area
#' spanned by 'dB' must fully contain all of the source support data, i.e., no
#' elements of 'spatialData' can lie outside of the boundary defined by 'dB'.
#' If 'spatialData' is multi-resolution 'dB' must be provided as the integer
#' element of 'spatialData' to be used as the finest resolution or as a 
#' SpatialPolygons object defining the finest resolution.
#' 
#' Input 'gbf' allows users to specify radial basis function beyond
#' the internally implemented bi-square, Wendland, and Gaussian radial functions. 
#' ('bisquare', 'wendland', 'gaussian')
#' If user provides
#' a function, the following formal arguments are required: 
#' crd - the coordinates at which the basis functions are to be evaluated; 
#' knots - the knots of the basis functions; and 
#' w - the scaling factor for the basis function. 
#' The function must return a matrix of dimension 
#' \{length(crd) x nrow(knots)\}.
#' 
#' For completeness, the bi-square function implemented in the package
#' is of the form
#' \deqn{\Psi_{j}(s) \equiv \left\{
#' \begin{array}{rl} 
#' \{1 - (||\mathrm{crd} - \mathrm{knots}_j||/w)^2\}^2 & \mathrm{if} ~ ||\mathrm{crd}-\mathrm{knots}_j|| \leq w ,\\
#' 0 & \mathrm{otherwise} \end{array} \right. .}
#' 
#' The Wendland function is
##' \deqn{ \Psi_{j}(s) \equiv \left\{
#' \begin{array}{rl} 
#'  \{ 1 - d_{j}(s)\}^6 \{35 d_{j}(s)^2 + 18 d_j(s) + 3\}/3 ~ 0 \leq d_{j} \leq 1, \\
#' 0 & \mathrm{otherwise} \end{array} \right. ,}
#' where 
#' \deqn{ d_{j}(s) = ||s - c_j||/w.}
#'
#' The Gaussian radial function is
#' \deqn{ \Psi_{j}(s) =
#'  \exp\{- \frac{1}{2} (||s - c_j||/w)^2\}.}
#'
#' For clarity, the default MI prior function returns
#' \deqn{K^{-1} = R_B^{-1} \mathcal{A} \{ Q^{'}_{B} (I-A)Q_{B}\}R_B^{-1},}
#' as defined in the supplemental section of the original manuscript.
#' 
#' This package implements methods of the \code{ff} package in
#' scenarios when memory size is of concern. The \code{ff} package
#' stores variables on the disk rather than in RAM. Though it is 
#' highly optimized, this choice does increase computation time.
#' Intermediate variables used only internally are stored in the
#' temp directory specified in input 'ffdir.' To trigger
#' this implementation, 'ffdir' must be set.
#' 
#' @name region
#' @rdname region
#' 
#' @param spatialData SpatialXXDataFrame as defined by 
#'   package sp or a list of said objects. Currently, this implementation 
#'   is limited to use of a
#'   SpatialPolygonsDataFrame or a SpatialPointsDataFrame.
#'   For multi-resolution, a list of said SpatialXXDataFrames. 
#'
#' @param response A numeric vector, character, or list of such. 
#'   If numeric, the value of interest; the vector must contain the response 
#'   for all spatialData. Specifically, if 'spatialData' is a list, vector must 
#'   include the value
#'   of interest for all data of element [[1]] followed by that for all 
#'   data of element [[2]], etc.  
#'   If a character, the column header of the data slot that holds the  
#'   value of interest; if a list is provided in 'spatialData' 
#'   each Spatial object must have the specified column header. If provided
#'   as a list, the elements must correspond to the elements of the 
#'   spatialData list.
#'
#' @param sigmavar A numeric vector/matrix, character, or list of said objects. 
#'   The survey variance.
#'   If numeric, the vector is the diagonal elements of the variance matrix
#'   and it must contain the variance 
#'   for all areal units of 'spatialData'. Specifically, if 'spatialData' is a 
#'   list, the vector must 
#'   include the variance for all data of element [[1]] followed by that for all 
#'   data of element [[2]], etc.  Similarly, if a matrix.
#'   If a character, the column header of the data slot that holds the  
#'   diagonal elements of the variance; if a list is provided in 'spatialData' 
#'   each Spatial object must have the specified column header. If provided
#'   as a list, the elements must correspond to the elements of the 
#'   'spatialData' list.
#'
#' @param dC A SpatialPolygons object or a vector defining the 
#'   desired clustering of the spatial data. If vector, it must be of the
#'   length of dB and contain the cluster id for each areal unit of dB.
#'
#' @param \dots ignored. Used to require named input for remaining formals.
#'
#' @param cage A logical. If TRUE, CAGE is estimated. If FALSE, DCAGE is 
#'    estimated.
#'
#' @param w A numeric. The scaling factor for radial basis functions. 
#'   See details for further information.
#'
#' @param longlat A logicial. TRUE indicates that spatialData is long/lat 
#'   coordinate system.
#'
#' @param nGibbs An integer. The total number of Gibbs samples generated.
#'
#' @param nBurn An integer. The first sample accepted in the 
#'   Gibbs sampling algorithm.
#'
#' @param nThin An integer. Keep every nThin sample of the Gibbs sampling
#'   algorithm once the burn specification has been satisfied.
#'
#' @param x NULL, character, or numeric. The covariates for the large scale
#'   variability model. If NULL, an intercept only model is assumed. 
#'   If numeric (matrix), the covariates; the object must contain the covariates 
#'   for all spatialData. Specifically, if 'spatialData' is a list, vector must 
#'   include the covariates for all data of element [[1]] followed by that for
#'   all data of element [[2]], etc.  
#'   If a character, the column header(s) of the data slot that holds the  
#'   covariates; if a list is provided in 'spatialData' 
#'   each Spatial object must have the specified column header. If provided
#'   as a list, the elements must correspond to the elements of the 
#'   spatialData list.
#'
#' @param sigma2_beta A numeric. The variance of the large scale variability
#'   model. 
#'
#' @param sigma2_xi A numeric. The variance of the fine scale variability. 
#'  Defaults to inverse Gamma.
#'
#' @param lambda A numeric. The r eigenvalues of the prior
#'   distribution on the \{r x r\} covariance matrix Q.
#'
#' @param dB NULL, integer, or a SpatialPolygons object defining the 
#'   finest resolution spatial object to be used for inference. If NULL,  
#'   'spatialData' is the finest resolution  
#'   considered. Note that if 'spatialData' is a list, 
#'   'dB' must be specified as either the element of 'spatialData' to be used
#'   as the finest resolution or as a SpatialPolygons object.
#'
#' @param gbf NULL, function or function name. The function to use to calculate
#'   the radial basis functions of the expansion of the   
#'   Obled-Creutin eigenfunction. The default bi-square function is available
#'   through this implementation. All others must be defined by user.  
#'   See details for further information. 
#'
#' @param Qprior A character. The Q prior. Must be one of {'MI' or 'wish'} 
#'   indicating the MI prior or the the wishart distribution respectively.
#'
#' @param nw An integer. The number of MC replicates to generate for estimating 
#'   the O-C eigenfunctions. If <=0, 'spatialData' must be or include 
#'   SpatialPoints data.
#'
#' @param knots A matrix or integer. If a matrix, the knots of the radial 
#'   basis functions. If an integer, the number of knots to generate using 
#'   fields::cover.design().
#'
#' @param dataScale A function. If the response of interest was scaled, the 
#'   function to undo the scaling. Used in plotting only.
#'
#' @param ffdir A character string. The directory in which ff object is to be 
#'   saved. See details for further information.
#'
#' @param nCore An integer. If using MC, the algorithm can be spread across
#'    nCore cores to expedite calculations. If nCore = 1L, no parallelization
#'    methods are used.
#'
#' @param wishartScale A numeric of NULL. If Qprior is "wishart", the value
#'   of the diagonal elements of the scale matrix provided to MCMCpack::riwish.
#'   Default value is 100. Input is ignored for all other values of Qprior.
#'
#' @param dMax A numeric. If finest resolution data is SpatialPoints and 
#'   Qprior is MI, the maximum distance at which points are considered
#'   adjacent. If not specified, it is taken to be upper boundary of the lowest
#'   decile.
#'
#' @param palette A character. The palette preference for plotting. The 
#'   palette is assumed to be viridis palette with possible values
#'   'viridis', 'magma', 'plasma', 'inferno', 'cividis'
#'
#' @param plot A logical. TRUE indicates that final plot will be generated.
#'
#' @param parallelLog A character object. A file name for logging
#'    parallel executions.
#'
#' @return A list.
#' \item{call     }{The original call structure}
#' \item{psi      }{A list containing the generating basis and OC weighting matrix.}
#' \item{CAGETrack}{The estimated CAGE/DCAGE for each dC cluster.}
#' \item{cluster  }{Clustering indices mapping dB to dC.}
#' \item{yOpt     }{The estimated value of 
#'                  interest at each Gibbs sample clustered according
#'                  to dC.}
#' \item{yFinest  }{The estimated value of 
#'                  interest at each Gibbs sample clustered according
#'                  to dB.}
#' \item{criterion}{"CAGE" or "DCAGE"}
#'
#' 
#' @references Bradley, J. R., Wikle, C. K., and Holan, S. H. (2017).
#'   Regionalization of Multiscale Spatial Processes using a Criterion
#'   for Spatial Aggregation Error. Journal of the Royal Statistical Society -
#'   Series B, 79, 815--832.
#' 
#' @examples
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
#' dC <- raster::rasterToPolygons(raster::raster(xmn = -1.25, xmx = 1.25,
#'                                               ymn = -1.25, ymx = 1.25,
#'                                               res = c(0.5,2.5),
#'                                               crs = "+proj=longlat +datum=WGS84"))
#'
#' knots <- cbind(c(-0.75, 0.0, 0.75, -0.75, 0.0, 0.75, -0.75, 0.0, 0.75),
#'                c(-0.75,-0.75, -0.75, 0.0, 0.0, 0.0, 0.75, 0.75, 0.75))
#'
#' res <- region(spatialData = dt,
#'               response = "x",
#'               sigmavar = rep(1, 25),
#'               dC = dC,
#'               nGibbs = 50L,
#'               nBurn = 10L,
#'               nThin = 1L,
#'               nw = 2000L,
#'               knots = knots)
#' 
#' 
#' 
#' @keywords cluster
#' @export
region <- function(spatialData,
                   response,
                   sigmavar,
                   dC,
                   ...,
                   cage = TRUE,
                   longlat = TRUE,
                   dB = NULL,
                   gbf = "bisquare",
                   w = NULL,
                   knots = 75L,
                   nw = 0L,
                   nGibbs = 10000L,
                   nBurn = 1000L,
                   nThin = 1L,
                   x = NULL,
                   sigma2_beta = 1e-4,
                   sigma2_xi = 1.0,
                   Qprior = "MI",
                   lambda = 1.0,
                   wishartScale = 100.0,
                   dMax = NULL,
                   ffdir = NULL,
                   nCore = 1L,
                   parallelLog = NULL,
                   plot = TRUE,
                   dataScale = NULL,
                   palette = "plasma") {

  Call <- match.call()

  # ensure that the number of cores is appropriately specified
  localCluster <- .check_nCore(nCore = nCore, parallelLog = parallelLog)

  ip <- .inputPrep(spatialData = spatialData,
                   response = response,
                   sigmavar = sigmavar,
                   rr = w,
                   longlat = longlat,
                   nGibbs = nGibbs,
                   nBurn = nBurn,
                   nThin = nThin,
                   x = x,
                   sigma2_beta = sigma2_beta,
                   sigma2_xi = sigma2_xi,
                   lambda = lambda,
                   dB = dB,
                   gbf = gbf,
                   Qprior = Qprior,
                   nw = nw,
                   knots = knots,
                   ffdir = ffdir,
                   localCluster = localCluster,
                   wishartScale = wishartScale,
                   dMax = dMax,
                   cage = cage)

  if (is.null(x = ip$dB)) {
    dB <- spatialData
  } else if (is.integer(x = ip$dB)) {
    dB <- spatialData[[ ip$dB ]]
  } else {
    dB <- ip$dB
  }

  # determine optimal CAGE
  optimal <- .dCCAGE(dB = dB,
                     yFinest = ip$yFinest,
                     ySource = ip$ySource,
                     finestOnSource = ip$finestOnSource,
                     sourceAreas = ip$sourceAreas,
                     finestAreas = ip$finestAreas,
                     dC = dC,
                     ffdir = ffdir)

  if (!is.null(localCluster)) {
    parallel::stopCluster(localCluster)
  }

  optimal[[ "psi" ]] <- ip$psi
  optimal[[ "yFinest" ]] <- ip$yFinest
  if (ip$criterion) optimal[[ "criterion" ]] <- "CAGE"
  if (!ip$criterion) optimal[[ "criterion" ]] <- "DCAGE"

  if (plot) {
    .mapIt(dB = dB,
           yFinest = optimal$yFinest,
           yOpt = optimal$yOpt,
           cage = optimal$CAGETrack,
           cluster = optimal$cluster$cluster,
           dataScale = dataScale,
           palette = palette,
           criterion = optimal$criterion)
  }


  if (is(object = optimal[[ "yFinest" ]], class2 = "ff_matrix")) {
    yFinest <- optimal$yFinest
    ffsave(list = "yFinest", file = file.path(ffdir, "yFinest"))
    optimal[[ "yFinest" ]] <- file.path(ffdir, "yFinest")
  }

  if (is(object = optimal[[ "yOpt" ]], class2 = "ff_matrix")) {
    yOpt <- optimal$yOpt
    ffsave(list = "yOpt", file = file.path(ffdir, "yOpt"))
  }

  optimal[[ "call" ]] <- Call

  class(optimal) <- "rcage"

  return( optimal )

}
