#
# Function is not exported
#
# Process input data to ensure appropriate structure and content and perform
# initial steps for method. Function verifies inputs, groups data to finest
# resolution, calculates knots, calls Qprior and generates Gibbs sample.
#
# @param \dots ignored. Used to require named input for remaining formals.
#
# @param spatialData A SpatialPointsDataFrame object, 
#   SpatialPolygonsDataFrame object or a list of said objects. The
#    source support data.
#
# @param response A numeric vector, character, or list of such. 
#   The outcome of interest.
#
# @param sigmavar A numeric vector, character, or list of such. 
#   The survey variance.
#
# @param rr A numeric object. The scaling factor for radial basis functions. 
#
# @param longlat A logical object. TRUE indicates that spatialData is in
#   long/lat coordinate system.
#
# @param nGibbs An integer object. The total number of Gibbs samples 
#   generated.
#
# @param nBurn An integer object. The first sample accepted in the 
#   Gibbs sampling algorithm.
#
# @param nThin An integer object. Keep every nThin sample of the Gibbs
#   sampling algorithm once the burn specification has been satisfied.
#
# @param x A character, numeric or NULL. The covariates for the large scale
#   variability model. If NULL, an intercept only model is assumed. 
#
# @param sigma2_beta A numeric object. The variance of the large scale
#   variability model. 
#
# @param sigma2_xi A numeric. The variance of the fine scale variability. 
#
# @param lambda A numeric vector. The r eigenvalues of the prior
#   distribution on the r x r covariance matrix Q.
#
# @param dB NULL or a SpatialPolygons object defining the 
#   finest resolution spatial object to be used for inference.
#
# @param gbf NULL, function or character. The function to use to calculate
#   the radial basis functions of the expansion of the   
#   Obled-Creutin eigenfunction.  
#
# @param Qprior A character object. The Q prior. Must be one of {'MI', 'wish'} 
#   indicating the MI prior or the Inverse Wishart distribution, respectively.
#
# @param nw An integer object or NULL. The number of MC replicates to generate
#   for estimating the components of the O-C basis. 
#
# @param knots A matrix or integer object. If a matrix, the x,y 
#   coordinates of the knots of the radial 
#   basis functions. If an integer, the number of knots to generate using 
#   fields::cover.design().
#
# @param ffdir A character or NULL. The directory in which ff objects are to
#   be saved.
#
# @param localCluster A cluster object or NULL.
#
# @param wishartScale A numeric or NULL. If Qprior is wish, the value
#   of the diagonal elements of the scale matrix provided to MCMCpack::riwish.
#   Default value is 100. Input is ignored for all other values of Qprior.
#
# @param dMax A numeric object. If finest resolution data is SpatialPoints and 
#   Qprior is MI, the maximum distance at which points are considered
#   adjacent. If not specified, it is taken to be the upper boundary of the
#   lowest decile.
#
# @param cage A logical object. If TRUE, CAGE is used to identify optimal.
#
# @return A list containing
#
# \item{gibbs         }{A list containing the parameters of the Gibbs sampling.}
# \item{yFinest       }{A matrix or ff_matrix of the modeled outcome on dB \{nB x nKept\}
# \item{finestAreas   }{A vector of the areas of the finest resolution areal units}
# \item{psi           }{A list of the components of the OC basis on dB}
# \item{dB            }{The finest resolution object or NULL}
# \item{gbfObj        }{The radial basis function object}
# \item{nw            }{The confirmed value for input nw}
# \item{ySource       }{A matrix or ff_matrix of the modeled outcome on the source}
# \item{finestOnSource}{Matrix connecting elements of finest resolution to source}
# \item{sourceAreas   }{A vector of the areas of the source support}
# \item{criterion     }{Logical indicating if CAGE or DCAGE}
# @name inputPrep
# @rdname inputPrep
#
#' @include checkSpatialData.R checkDB.R checkNW.R
#' @include verifyNumericVector.R verifyCovariates.R
#' @include verifyffdir.R GBFObj.R verifyQPrior.R
#' @include QPriorObj.R QPriorObj_MI.R QPriorObj_Wishart.R hMatrix.R
#' @include obledCruetinBasis.R gibbs.R generatingBasis.R
#
.inputPrep <- function(...,
                       spatialData,
                       response,
                       sigmavar,
                       rr,
                       longlat,
                       nGibbs,
                       nBurn,
                       nThin,
                       x,
                       sigma2_beta,
                       sigma2_xi,
                       lambda,
                       dB,
                       gbf,
                       Qprior,
                       nw,
                       knots,
                       ffdir,
                       localCluster,
                       wishartScale,
                       dMax,
                       cage) {

  # ensures that spatialData is a SpatialPointsDataFrame, 
  # SpatialPolygonsDataFrame or a list of said objects
  spatialData <- .checkSpatialData(spatialData = spatialData)

  # ensures appropriate nw value
  nw <- .checkNW(spatialData = spatialData, nw = nw)

  # ensures that dB is one of NULL, integer, or SpatialPolygons and 
  # that cage is set appropriately
  dB <- .checkDB(dB = dB, spatialData = spatialData, cage = cage)

  # cage TRUE indicates that SpatialPoints data is present and CAGE was
  # requested to be used to obtain optimal
  cage <- dB$cage

  # if dB provided as SpatialPolygons, returned 
  # if dB NULL at input and spatialData is a list, error
  # if dB NULL at input and spatialData is SpatialPolygons or SpatialPoints, NULL
  # if dB numeric and spatialData is a list, returns element id as provided
  dB <- dB$dB

  # extracts response from spatialData object or ensures
  # correct length of provided vector returned object is a
  # vector of length equal to the total number of areal units
  response <- .verifyNumericVector(spatialData = spatialData, 
                                   numVec = response,  
                                   nm = 'response')

  # extracts covariates from spatialData object, generates intercept only
  # model if NULL, or ensures correct length of provided vector/matrix
  # returned object is a matrix is dim 1 equal to the total number of areal
  # units
  x <- .verifyCovariates(spatialData = spatialData, x = x)

  # extracts variance from spatialData object or ensures
  # correct length of provided vector returned object is a
  # vector of length equal to the total number of areal units
  sigmavar <- .verifyNumericVector(spatialData = spatialData, 
                                   numVec = sigmavar,  
                                   nm = 'sigmavar')

  # ensures Gibbs sampler inputs are reasonable
  nBurn <- as.integer(x = round(x = nBurn, digits = 0L))
  if (nBurn <= 0L) {
    message("nBurn reset to 100")
    nBurn <- 100L
  }
  nThin <- as.integer(x = round(x = nThin, digits = 0L))
  if (nThin <= 0L) {
    message("nThin reset to 1")
    nThin <- 1L
  }

  # vector of indices of the gibbs samples to be kept
  gibbsKeep <- tryCatch(expr = seq(from = {nBurn+1L}, to = nGibbs, by = nThin),
                        condition = function(e){ stop(e$message) })
  message("Gibbs sampler will keep ", length(x = gibbsKeep), " samples")

  # ensures that there is *likely* sufficient memory to hold results
  # or that, if specified, the directory of storing objects exists
  # no object is returned
  tst <- .verifyffdir(ffdir = ffdir, 
                      nKept = length(x = gibbsKeep),
                      nR = length(x = response))

  # {nB x nSource}
  # relates the elements of the source data to the elements of dB
  # SpatialPoints data have 1's indicating in which element of dB they are
  # located. If the point lies on a boundary, it is randomly assigned to
  # one of the neighboring polygons.
  # SpatialPolygons data have fractions indicating how much of their
  # area is covered by each element of dB
  # this matrix is used to correlate the appropriate xi element(s) with the
  # source data
  finestOnSource <- .hMatrix(spatialData = spatialData, dB = dB)

  message("h matrix generated")

  # object of class GBFObj containing the weight, knots, and
  # function name for calculating the basis functions.
  # initialize this object
  gbfObj <- .verifyGBF(gbf = gbf, 
                       weight = rr, 
                       knots = knots, 
                       longlat = longlat,
                       spatialData = spatialData,
                       dB = dB)

  # ensures an appropriate QPriorObj indicating the function to be used to 
  # calculate Q
  Qprior <- .verifyQPrior(qPrior = Qprior)

  # generate components of OC basis functions
  # {nSource x r}
  message("Source Support")
  basis <- .obledCruetinBasis(spatialData = spatialData, 
                              dB = dB,
                              gbfObj = gbfObj, 
                              nw = nw, 
                              localCluster = localCluster,
                              verify = TRUE)

  # bandwidth of basis function might have been changed.
  gbfObj <- basis$gbfObj

  message("Finest Support")
  idB <- NA
  # generate basis on dB
  if (is.null(x = dB)) {
    # if dB is NULL using source support data is the finest resolution
    message("\tusing source support basis")
    basisdB <- basis$phiOC
  } else if (is.integer(x = dB)) {
    # if integer using an element of multi-resolution source support
    # generate radial basis functions on finest support (psi(A))
    basisdB <- .generatingBasis(spatialData = spatialData[[ dB ]],
                                nw = nw,
                                gbfObj = gbfObj,
                                db = min(10000L, nw),
                                localCluster = localCluster)

    # normalize
    basisdB <- basisdB %*% basis$OCnorm

    ortho <- t(x = basisdB) %*% basisdB
    for (i in 1L:ncol(x = basisdB)) {
       basisdB[,i] <- basisdB[,i] / sqrt(ortho[i,i])
    }
  } else {
    # generate radial basis functions on finest support (psi(A))
    basisdB <- .generatingBasis(spatialData = dB,
                                nw = nw,
                                gbfObj = gbfObj,
                                db = min(10000L, nw),
                                localCluster = localCluster)

    # normalize
    basisdB <- basisdB %*% basis$OCnorm

    ortho <- t(x = basisdB) %*% basisdB
    for (i in 1L:ncol(x = basisdB)) {
       basisdB[,i] <- basisdB[,i] / sqrt(ortho[i,i])
    }
  }

  # calculate Q inverse
  Qprior <- .qInv(qObj = Qprior,
                  basisdB = basisdB,
                  lambda = lambda,
                  scale = wishartScale,
                  spatialData = spatialData,
                  dB = dB,
                  dMax = dMax)

  # there may be circumstances where covariate information is available
  # but not outcome. Remove these cases from the Gibbs step
  naResponse <- is.na(x = response) | is.nan(x = response)

  if (sum(naResponse) > 0L) {
    message("excluding ", sum(naResponse), 
            " cases from Gibbs due to incomplete response data")
  }

  notna <- !naResponse

  #{nSource' x nB}
  spatialOnFinest <- Matrix::Matrix(data = t(x = finestOnSource[,notna,drop=FALSE]), 
                                    sparse = TRUE)

  # eliminate empty cases of dB
  isEmpty <- colSums(as.matrix(x = spatialOnFinest) > 1e-8) == 0L

  if (any(isEmpty)) {
    message("excluding ", sum(isEmpty), 
            " elements of D_B from Gibbs step for lack of data")
    #{nSource' x nB'}
    spatialOnFinest <- spatialOnFinest[,!isEmpty]
  }

  # Gibbs sampling for beta, eta, and xi
  # returns a list object with 
  #  beta { p x nGibbs } 
  #  eta { r x nGibbs }
  #  xi { nB' x nGibbs } 
  if (is.matrix(x = sigmavar)) {
    tsigmavar <- sigmavar[notna,notna,drop=FALSE]
  } else {
    tsigmavar <- sigmavar[notna]
  }

  gibbsSamples <- .gibbs(Z = response[notna], #{nSource'}
                         X = x[notna,,drop=FALSE], #{nSource' x nBeta}
                         H = spatialOnFinest, #{nSource' x nB}
                         psi = basis$phiOC[notna,,drop=FALSE], #{nSource' x r}
                         sigma2_eps = tsigmavar,#{nSource'}
                         sigma2_beta = sigma2_beta, #{np}
                         sigma2_xi = sigma2_xi, #{1}
                         qObj = Qprior,
                         gibbsKeep = gibbsKeep,
                         ffdir = ffdir) #{nB}

  message("finished Gibbs sampling")

  # add back empty cells of dB to xi
  # {nB x nGibbs}
  if (any(isEmpty)) {
    xi <- .makeStorage(ffdir = ffdir, 
                       nrow = nrow(x = finestOnSource), 
                       ncol = ncol(x = gibbsSamples$xi))

    for (i in 1L:ncol(x = gibbsSamples$xi)) {
      xi[!isEmpty,i] <- gibbsSamples$xi[,i]
    }

    gibbsSamples$xi <- xi
  }

  ## yFinest {nB x nGibbs}
  yFinest <- .yFinest(ffdir = ffdir, 
                      gibbsSamples = gibbsSamples, 
                      x = x,
                      spatialData = spatialData,
                      dB = dB,
                      basisdB = basisdB,
                      fos = finestOnSource)

  if (is.null(x = dB)) {

    # if dB is NULL using source support data as finest resolution

    finestAreas <- .getArea(spatialData = spatialData, byid = TRUE)

  } else if (is.numeric(x = dB)) {

    # if dB is an integer, using an element of the source support as dB

    finestAreas <- .getArea(spatialData = spatialData[[ dB ]], byid = TRUE)

  } else {

    finestAreas <- .getArea(spatialData = dB, byid = TRUE)

  }

  ySource <- .yCageSource(ffdir = ffdir, 
                          gibbsSamples = gibbsSamples, 
                          x = x,
                          spatialData = spatialData,
                          cage = cage,
                          basis = basis,
                          H = t(x = finestOnSource))

  # indicator identifies the components of the source data used to 
  # calculate cage/DCAGE
  # {nSourceCage x nB}
  fos <- ySource$indicator %*% t(x = finestOnSource)

  # need to re-evaluate isEmpty to include elements removed for being na
  # {nB}
  isEmpty <- rowSums(as.matrix(x = finestOnSource) > 1e-8) == 0L

  return( list("gibbs" = gibbsSamples,
               "yFinest" = yFinest, # {nB x nGibbs}
               "finestAreas" = finestAreas, # {nB}
               "psi" = basis,
               "dB" = dB,
               "gbfObj" = gbfObj,
               "nw" = nw,
               "ySource" = ySource$ySource, # {nGibbs x nSource'}
               "finestOnSource" = fos, # {nSource' x nB}
               "sourceAreas" = ySource$areas, # {nSource'}
               "criterion" = cage,
               "isEmpty" = isEmpty) ) #{nB}

}
