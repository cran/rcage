#  Determine best clustering based on minimum CAGE/DCAGE
#
#  Function not exported
#
#  @param ... ignored.
#
#  @param dB A SpatialPolygons object or NULL. The finest resolution {nB}
#
#  @param yFinest A matrix object. The modeled outcome from Gibbs sampling
#    on finest resolution.
#    {nB x nKept}
#
#  @param ySource A matrix object. The modeled outcome from Gibbs sampling
#    on source support.
#    {nKept x nSource}
#
#  @param finestOnSource A matrix object. Connects source support to
#    finest resolution.
#    {nSource x nB}
#
#  @param sourceAreas A vector object. The area of each source areal unit.
#
#  @param finestAreas A vector object. The area of each finest resolution 
#    areal unit.
#
#  @param gL A integer object. The smallest number of areal units to 
#    consider.
#
#  @param gU A integer object. The largest number of areal units to 
#    consider. 
#
#  @param cMethod A character object. The clustering method.
#
#  @param alpha A numeric object. Between 0 and 1; mixing parameter gives
#    the relative importance of D0 compared to D1.
#
#  @param localCluster A cluster object or NULL. A cluster for parallel.
#
#' @include getArea.R minCAGE.R minCAGESH.R yOpt.R
#' @importFrom rgeos gCentroid
.optimalDCAGE <- function(...,
                          dB,
                          yFinest,
                          ySource,
                          finestOnSource,
                          sourceAreas,
                          finestAreas,
                          gL,
                          gU,
                          cMethod,
                          alpha,
                          localCluster) {

  message("determining clustering that minimizes CAGE/DCAGE")

  # centroids of areal units
  # if dB is spatialPoints, this returns the coordinates of the points
  centroids <- tryCatch(expr = sp::coordinates(obj = rgeos::gCentroid(spgeom = dB, 
                                                                      byid = TRUE)),
                        error = function(e) {
                                stop("unable to identify centroids of finest ",
                                     "resolution spatial data\n", e$message, 
                                     call. = FALSE)
                              })

  if (cMethod == 'kmeans') {

    minResult <- .minCAGE(centroids = centroids,
                          yFinest = yFinest,
                          ySource = ySource,
                          finestOnSource = finestOnSource,
                          sourceAreas = sourceAreas,
                          finestAreas = finestAreas,
                          gL = gL,
                          gU = gU,
                          localCluster = localCluster)

  } else {

    minResult <- .minCAGESH(dB = dB,
                            centroids = centroids,
                            yFinest = yFinest,
                            ySource = ySource,
                            finestOnSource = finestOnSource,
                            sourceAreas = sourceAreas,
                            finestAreas = finestAreas,
                            gL = gL,
                            gU = gU,
                            alpha = alpha,
                            localCluster = localCluster)
  }

  # area weighted value of interest in each cluster for Gibbs samples
  YOpt <- .yOpt(idxit = minResult$cluster$cluster, 
                y = yFinest,
                areas = finestAreas)

  minResult[[ "yOpt" ]] <- YOpt

  return( minResult )

}
