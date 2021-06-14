# @param ... Ignored.
#
# @param dB A SpatialPolygons object, integer, or NULL. The finest resolution.
#
# @param yFinest A matrix or ff_matrix object. The modeled outcome of 
#   interest on the finest resolution. {nB x nGibbs}
#
# @param ySource A matrix object. The modeled outcome from Gibbs sampling
#   on source support.
#   {nGibbs x nSource}
#
# @param finestOnSource A matrix object. Connects source support to
#   finest resolution.
#   {nSource x nB}
#
# @param sourceAreas A vector object. The area of each areal unit in source.
#
# @param finestAreas A vector object. The area of each areal unit in finest.
#
# @param dC A SpatialPolygons object or vector. The desired clustering.
#
# @param ffdir A character or NULL. ff director
#
# @return A list containing
#   \item{CAGETrack}{The estimated CAGE/DCAGE for each dC cluster.}
#   \item{cluster  }{Clustering indices mapping dB to dC.}
#   \item{yOpt     }{The estimated value of 
#                    interest at each Gibbs sample clustered according
#                    to dC.}
#
#' @include hMatrix.R cage.R yOpt.R
.dCCAGE <- function(..., 
                    dB,
                    yFinest,
                    ySource,
                    finestOnSource,
                    sourceAreas,
                    finestAreas,
                    dC,
                    ffdir) {

  message("calculating CAGE/DCAGE on dC")

  # indicator of overlap of dB and dC {nCluster x nB}
  if (is(object = dC, class2 = "Spatial")) {

    hm <- .hMatrix(spatialData = dB, dB = dC, normalize = FALSE)

    if (!any(hm > 0)) stop("dC does not overlay dB", call. = FALSE)

    # cluster id for each areal unit in dB {nB}
    dCIDs <- unname(obj = drop(x = t(x = hm > 0.0) %*% 1L:length(x = dC)))

  } else {

    if (length(x = dC) != length(x = dB)) stop("dim error in dC", call. = FALSE)

    dCIDs <- dC
  }


  i1 <- NULL
  i2 <- NULL

  # DCAGE on dC 
  cageForIndex <- .cage(yFinest = yFinest,
                        ySource = ySource,
                        finestOnSource = finestOnSource,
                        idxit = dCIDs,
                        sourceAreas = sourceAreas,
                        finestAreas = finestAreas)

  # area weighted response in each cluster for Gibbs samples
  # {nCluster x nGibbs}
  yOpt <- .yOpt(idxit = dCIDs, 
                y = yFinest,
                areas = finestAreas)

  return( list("CAGETrack" = cageForIndex,
               "cluster"   = list("cluster" = dCIDs),
               "yOpt"      = yOpt) )

}
