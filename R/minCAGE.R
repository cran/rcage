#  Determine Best Clustering Based on Minimum CAGE
#
#  function is not exported
#
#  @param centroids A "vector" object. The centroids. {nB}
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
#  @param gL An integer object. The smallest number of areal units to 
#    consider.
#
#  @param gU An integer object. The largest number of areal units to 
#    consider.
#
#  @param localCluster A cluster object or NULL. The cluster if parallel.
#
#' @importFrom stats kmeans
#' @importFrom parallel parApply parLapply
#' @include cage.R
.minCAGE <- function(...,
                     centroids,
                     yFinest,
                     ySource,
                     finestOnSource,
                     sourceAreas,
                     finestAreas,
                     gL,
                     gU,
                     localCluster) {

  if (!is.null(x = localCluster)) {

    message("obtaining kmeans in parallel")

    result <- parallel::parSapply(cl = localCluster, 
                                  X = 1L:ncol(x = yFinest),
                                  FUN = .allClusterKmeans,
                                  centroids = centroids,  
                                  gL = gL,  
                                  gU = gU,
                                  yFinest = yFinest,
                                  ySource = ySource,
                                  finestOnSource = finestOnSource,
                                  sourceAreas = sourceAreas,
                                  finestAreas = finestAreas)

    tst <- which(result <= {min(result)+1e-8}, arr.ind = TRUE)

    whichGibbs <- tst[1L,2L]
    whichG <- tst[1L,1L]

    minCAGE <- result[whichG, whichGibbs]

  } else {

    minCAGE <- Inf
    whichGibbs <- 0L
    whichG <- 0L

    message("kmeans for Gibbs sample:")

    for (i in 1L:ncol(x = yFinest)) {

      message(i, " ", appendLF = FALSE)
      if (i %% 10L == 0L) message("")

      result <- .allClusterKmeans(x = i,  
                                  centroids = centroids,  
                                  gL = gL,  
                                  gU = gU,
                                  yFinest = yFinest,
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

  if (is.infinite(x = minCAGE)) stop("infinite result in CAGE/DCAGE", 
                                     call. = FALSE)

  Jmat <- cbind(centroids, yFinest[,whichGibbs])

  IDX <- tryCatch(expr = stats::kmeans(x = Jmat,
                                       centers = {gL:gU}[whichG],
                                       iter.max = 10000L),
                  error = function(e) {
                            stop("unable to obtain clustering\n",
                                 e$message, call. = FALSE)
                          })

  cageForIndex <- .cage(yFinest = yFinest,
                        ySource = ySource,
                        finestOnSource = finestOnSource,
                        idxit = IDX$cluster,
                        sourceAreas = sourceAreas,
                        finestAreas = finestAreas)

  return( list("minCAGE" = mean(x = cageForIndex),
               "CAGETrack" = cageForIndex,
               "cluster" = IDX) )
}

#  function to obtain DCAGE/CAGE at all clusterings for a single Gibbs samples
#
#  @param x A vector object. The Gibbs sample to consider.
#
#  @param ... ignored
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
#    {nS x nB} or {nA x nB}
#
#  @param sourceAreas A vector object. The areas of areal units.
#
#  @param centroids A matrix object. The centroids of the areal units.
#
#  @param gL An integer object. The smallest number of clusters.
#
#  @param gU An integer object . The largest number of clusters.
#
#' @importFrom stats kmeans
.allClusterKmeans <- function(x,
                              ...,
                              yFinest,
                              ySource,
                              finestOnSource,
                              sourceAreas,
                              finestAreas,
                              centroids,  
                              gL,  
                              gU) {

  # Jmatrix for use in kmeans
  # { n x 3 }

  Jmat <- cbind(centroids, yFinest[,x])

  IDX <- tryCatch(expr = stats::kmeans(x = Jmat,
                                       centers = gU,
                                       iter.max = 10000L),
                  error = function(e) {
                            stop("unable to obtain clustering\n",
                                 e$message, call. = FALSE)
                          })

  cageForIndex <- .cage(yFinest = yFinest,
                        ySource = ySource,
                        finestOnSource = finestOnSource,
                        sourceAreas = sourceAreas,
                        finestAreas = finestAreas,
                        idxit = IDX$cluster)

  cage <- mean(x = cageForIndex)

  i <- gU - 1L

  while (i >= gL) {

    # cluster data into i clusters
    IDX <- tryCatch(expr = stats::kmeans(x = Jmat,
                                         centers = IDX$centers[1L:i,1L:3L,drop=FALSE],
                                         iter.max = 10000L),
                    error = function(e) {
                              msg <- paste0("unable to obtain clustering\n",
                                            e$message)
                              stop(msg)
                            })

    cageForIndex <- .cage(yFinest = yFinest,
                          ySource = ySource,
                          finestOnSource = finestOnSource,
                          sourceAreas = sourceAreas,
                          finestAreas = finestAreas,
                          idxit = IDX$cluster)

    # average gamma across all clusters
    cage <- c(mean(x = cageForIndex), cage)

    i <- i - 1L

  }

  return( cage )
}
