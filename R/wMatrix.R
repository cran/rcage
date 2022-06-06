# spatialData	nw	dB	W calculated using
# pts		NULL	--	spatialData
# pts		nMC	NULL	nMC samples from convex hull of pts
# pts		nMC	poly	nMC samples from dB
# poly		NULL	NULL	nMC samples from spatialData (this should not happen
#                               if checkNW is working properly)
# poly		nMC	--	nMC samples from spatialData
# list		NULL	--	spatialPoints element
# list		nMC	poly	nMC sample of dB
# list		nMC	int	nMC samples of dB element

#' @importFrom sp coordinates
#' @include getArea.R
.coreFunc_wMatrix <- function(..., pts, gbfObj, spatialData) {

  # generate basis on samples points {nw x r}
  basis <- .basis(crd = sp::coordinates(obj = pts), gbfObj = gbfObj)

  # area of the full spatial domain {1}
  area <- .getArea(spatialData = spatialData, byid = FALSE)

  W <- area * crossprod(x = basis) / nrow(x = basis)

  # {r x r}
  # A m^{-1} sum_(k=1)^{m} g_i(s_k) g_j(s_k)
  return( W )

}

#  W Matrix
#
#  method is not exported
#
#  @param spatialData A SpatialPoints object, SpatialPolygons object, or 
#   list of said objects. The source support.
#
#  @param dB A SpatialPolygons object, NULL, or numeric.  
#   The finest resolution.
#
#  @param nw An integer object or NULL. The number of MC replicates.
#
#  @param ... ignored.
#
#  @param gbfObj A GBFObj object. Information about basis function.
#
#  @param verbose A logical object. If TRUE provide messages.
#
setGeneric(name = ".wMatrix",
           def = function(spatialData, dB, nw, ...) { 
                   standardGeneric(".wMatrix") 
                 })

# default method results in error
setMethod(f = ".wMatrix",
          signature = c(spatialData = "ANY",
                        dB = "ANY",
                        nw = "ANY"),
          definition = function(spatialData, dB, nw, ...) { 
                         stop("not allowed") 
                       })

#  SpatialPoints are sufficient to estimate Psi
#' @importFrom sp coordinates
setMethod(f = ".wMatrix",
          signature = c(spatialData = "SpatialPoints",
                        dB = "ANY",
                        nw = "NULL"),
          definition = function(spatialData, dB, nw, ..., gbfObj, verbose = TRUE) {

              if (verbose) message("\tcalculating W using SpatialPoints data")

              return( .coreFunc_wMatrix(pts = sp::coordinates(obj = spatialData), 
                                        gbfObj = gbfObj,  
                                        spatialData = spatialData) )

            })


#  SpatialPoints are provided but are not sufficient to estimate W
#  User did not provide finest resolution object; use convex hull
#' @importFrom rgeos gConvexHull
#' @importFrom sp spsample
setMethod(f = ".wMatrix",
          signature = c(spatialData = "SpatialPoints",
                        dB = "NULL",
                        nw = "numeric"),
          definition = function(spatialData, dB, nw, ..., gbfObj, verbose = TRUE) {

              if (verbose) {
                message("\tconvex hull of spatialData used as sample region",
                        " for MC estimate of W")
              }

              dB <- rgeos::gConvexHull(spgeom = spatialData)

              # sample spatial data {nw x 2}
              mcPoints <- sp::spsample(x = dB, n = nw, type = "random")

              return( .coreFunc_wMatrix(pts = mcPoints, 
                                        gbfObj = gbfObj,  
                                        spatialData = dB) )

            })

#  SpatialPoints are provided but are not sufficient to estimate W
#  User provided finest resolution object; sample from dB
#' @importFrom sp spsample
setMethod(f = ".wMatrix",
          signature = c(spatialData = "SpatialPoints",
                        dB = "SpatialPolygons",
                        nw = "numeric"),
          definition = function(spatialData, dB, nw, ..., gbfObj, verbose = TRUE) {

              if (verbose) {
                message("\tdB used as sample region for MC estimate of W")
              }

              # sample spatial data {nw x 2}
              mcPoints <- sp::spsample(x = dB, n = nw, type = "random")

              return( .coreFunc_wMatrix(pts = mcPoints, 
                                        gbfObj = gbfObj,  
                                        spatialData = dB) )

            })

#  MC estimate of W when SpatialPoints are not provided
#  Regardless of finest resolution specification, sample polygons data
#' @importFrom sp spsample
setMethod(f = ".wMatrix",
          signature = c(spatialData = "SpatialPolygons",
                        dB = "ANY",
                        nw = "numeric"),
          definition = function(spatialData, dB, nw, ..., gbfObj, verbose = TRUE) {

              if (verbose) {
                message("\tspatialData used as sample region for MC estimate of W")
              }

              # sample spatial data {nw x 2}
              mcPoints <- sp::spsample(x = spatialData,  
                                       n = nw, 
                                       type = "random")

              return( .coreFunc_wMatrix(pts = mcPoints, 
                                        gbfObj = gbfObj,  
                                        spatialData = spatialData) )

            })

# if a list containing multi-resolution data is provided 
#   and nw is NULL, SpatialPoints data is sufficient
#   to estimate W or psi(A). 
setMethod(f = ".wMatrix",
          signature = c(spatialData = "list",
                        dB = "ANY",
                        nw = "NULL"),
          definition = function(spatialData, dB, nw, ...) {

              for (i in 1L:length(x = list)) {
                if (is(object = spatialData[[ i ]],
                       class2 = "SpatialPoints")) {
                  return( .wMatrix(spatialData = spatialData[[ i ]],
                                   dB = NULL,
                                   nw = NULL, ...) )
                }
              }

            })

# if a list containing multi-resolution data is provided 
#   and nw is numeric, SpatialPoints data are not sufficient
#   to estimate W or psi(A). Use dB as the sample region; dB provided
#   by user.
#' @importFrom sp spsample
setMethod(f = ".wMatrix",
          signature = c(spatialData = "list",
                        dB = "SpatialPolygons",
                        nw = "numeric"),
          definition = function(spatialData, dB, nw, ..., gbfObj, verbose = TRUE) {

              if (verbose) {
                message("\tdB used as sample region for MC estimate of W")
              }

              # sample spatial data {nw x 2}
              mcPoints <- sp::spsample(x = dB, n = nw, type = "random")

              return( .coreFunc_wMatrix(pts = mcPoints, 
                                        gbfObj = gbfObj,  
                                        spatialData = dB) )

            })

# if a list containing multi-resolution data is provided 
#   and nw is numeric, SpatialPoints data are not sufficient
#   to estimate W or psi(A). Use dB as the sample region; dB taken from
#   source support.
#' @importFrom sp spsample
setMethod(f = ".wMatrix",
          signature = c(spatialData = "list",
                        dB = "numeric",
                        nw = "numeric"),
          definition = function(spatialData, dB, nw, ..., gbfObj, verbose = TRUE) {

              if (verbose) {
                message("\tdB used as sample region for MC estimate of W")
              }

              # sample spatial data {nw x 2}
              mcPoints <- sp::spsample(x = spatialData[[ dB ]],  
                                       n = nw, 
                                       type = "random")

              return( .coreFunc_wMatrix(pts = mcPoints, 
                                        gbfObj = gbfObj,  
                                        spatialData = spatialData[[ dB ]]) )

            })
