#  Determine Best Clustering Based on Minimum CAGE/DCAGE Using Hierarchical
#
#  @param dB A SpatialPolygons object or NULL. The finest resolution.
#
#  @param centroids A vector object. The centroids. {nB}
#
#  @param yFinest A ff_matrix or matrix object. The outcome of interest
#    on the finest resolution. 
#    {nB x nKept}
#
#  @param ySource A matrix object. The outcome of interest on the source
#    support. 
#    {nKept x nS} or {nKept x nA}
#
#  @param finestOnSource A matrix object. Membership matrix linking
#    source support to finest resolution. 
#    {nA x nB} or {nS x nB}
#
#  @param sourceAreas A vector object. Areas of source areal units.
#    {nB}
#
#  @param gL An integer object. The smallest number of areal units to 
#    consider.
#
#  @param gU An integer object. The largest number of areal units to 
#    consider.
#
#  @param localCluster A cluster object or NULL. The cluster if parallel.
#
#  @param alpha A numeric object. Between 0 and 1; mixing parameter gives
#    the relative importance of D0 compared to D1.
#
#' @importFrom rgeos gCentroid gTouches
#' @importFrom stats as.dist dist cutree
#' @importFrom ClustGeo hclustgeo
#' @importFrom parallel parApply
.minCAGESH <- function(...,
                       dB,
                       centroids,
                       yFinest,
                       ySource,
                       finestOnSource,
                       sourceAreas,
                       finestAreas,
                       gL,
                       gU,
                       localCluster,
                       alpha) {

  # adjacency matrix
  D1 <- sp::spDists(x = dB)
  D1 <- D1 / max(D1)
  D1 <- stats::as.dist(m = D1)

  minCAGE <- Inf

  if (!is.null(x = localCluster)) {

    message("obtaining hclustgeo in parallel")

    result <- parallel::parSapply(cl = localCluster, 
                                  X = 1L:ncol(x = yFinest), 
                                  MARGIN = 2L,
                                  FUN = .allClusterSH,
                                  yFinest = yFinest,
                                  centroids = centroids,  
                                  gL = gL,  
                                  gU = gU,  
                                  D1 = D1,
                                  alpha = alpha,
                                  ySource = ySource,
                                  finestOnSource = finestOnSource,
                                  sourceAreas = sourceAreas,
                                  finestAreas = finestAreas)

    tst <- which(result <= {min(result)+1e-8}, arr.ind = TRUE)

    whichGibbs <- tst[1L,2L]
    whichG <- tst[1L,1L]

    minCAGE <- result[whichG, whichGibbs]

  } else {

    message("hclustgeo for Gibbs sample:")

    minCAGE <- Inf
    whichGibbs <- 0L
    whichG <- 0L

    for (i in 1L:ncol(x = yFinest)) {

      message(i, " ", appendLF = FALSE)
      if (i %% 10L == 0L) message("")

      result <- .allClusterSH(x = i,  
                              yFinest = yFinest,
                              centroids = centroids,  
                              gL = gL,  
                              gU = gU,  
                              D1 = D1,
                              alpha = alpha,
                              ySource = ySource,
                              finestOnSource = finestOnSource,
                              sourceAreas = sourceAreas, 
                              finestAreas = finestAreas)

      if (min(result) < minCAGE) {
        whichGibbs <- i
        whichG <- which.min(x = result)
        minCAGE <- result[whichG]
      }
    }

    message("")
  }

  

  if (is.infinite(x = minCAGE)) stop("Inf result in CAGE/DCAGE")

  # Jmatrix for use in kmeans
  # { n x 2 }
  Jmat <- cbind(centroids, yFinest[,whichGibbs])

  # cluster data
  IDX <- tryCatch(expr = ClustGeo::hclustgeo(D0 = dist(x = Jmat),
                                             D1 = D1,
                                             alpha = alpha),
                  error = function(e) {
                            msg <- paste0('unable to obtain clustering\n', e$message)
                            stop(msg)
                          })

  # cut tree down to specified level
  cluster <- cutree(tree = IDX, k = {gL:gU}[whichG])

  # calculate CAGE/DCAGE
  cageForIndex <- .cage(yFinest = yFinest,
                        ySource = ySource,
                        finestOnSource = finestOnSource,
                        sourceAreas = sourceAreas,
                        finestAreas = finestAreas,
                        idxit = cluster)

  return( list("minCAGE" = mean(x = cageForIndex),
               "CAGETrack" = cageForIndex,
               "cluster" = list("hclust" = IDX,
                                "cluster" = cluster)) )
}

# function to obtain DCAGE/CAGE at all clusterings for a single
# Gibbs sample.
# Input j is the index for the Gibbs sample to consider
.allClusterSH <- function(x, 
                          yFinest, 
                          ySource,
                          finestOnSource, 
                          sourceAreas,
                          finestAreas,
                          centroids,  
                          gL,  
                          gU,  
                          inCluster, 
                          D1,   
                          alpha, ...) {

  # Jmatrix for use in kmeans
  # { n x 2 }
  Jmat <- cbind(centroids, yFinest[,x])

  # cluster data
  IDX <- tryCatch(expr = ClustGeo::hclustgeo(D0 = dist(x = Jmat),
                                             D1 = D1,
                                             alpha = alpha),
                  error = function(e) {
                            stop('unable to obtain clustering\n', 
                                 e$message, call. = FALSE)
                          })

  cage <- NULL

  for (i in gL:gU) {

    # cut tree down to specified level
    cluster <- cutree(tree = IDX, k = i)

    # calculate CAGE/DCAGE
    cageForIndex <- .cage(yFinest = yFinest,
                          ySource = ySource,
                          finestOnSource = finestOnSource,
                          sourceAreas = sourceAreas,
                          finestAreas = finestAreas,
                          idxit = cluster)

    # average gamma across all clusters
    newCAGE <- mean(x = cageForIndex)

    cage <- c(cage, newCAGE)

  }

  return( cage )
}
