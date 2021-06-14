#  Calculate CAGE or DCAGE
#
#  function is not exported
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
# @param finestAreas A vector object. The area of each areal unit in finest.
#    {nB}
#
# @param sourceAreas A vector object. The area of each areal unit in source.
#    {nSource}
#
# @param idxit A vector object.  The cluster to which each aggregated
#   spatial unit in finest resolution belongs. {nB}
#
#' @importFrom stats aggregate
#' @import Matrix
.cage <- function(..., yFinest, 
                       ySource, 
                       finestOnSource, 
                       sourceAreas, 
                       finestAreas,
                       idxit) {

  # unique cluster identifiers 
  # {nC}
  unqClusters <- sort(x = unique(x = idxit))

  # number of unique clusters
  nC <- length(x = unqClusters)

  # number of areal units in finest resolution
  nB <- length(x = idxit)

  # cluster to which each areal unit is a member {nB}
  idx <- match(x = idxit, table = unqClusters)

  # matrix linking element of dB to each cluster {nC x nB}
  clusterOnDB <- Matrix::Matrix(data = 0.0, nrow = nC, ncol = nB, sparse = TRUE)
  clusterOnDB[cbind(idx,1L:nB)] <- 1L

  # total area for each cluster {nCluster} 
  # sums the areas of the finest resolution units in the cluster
  totalAreas <- drop(x = as.matrix(x = clusterOnDB %*% finestAreas))

  # convert indicator matrix to proportion
  clusterOnDB2 <- clusterOnDB / totalAreas

  #{nC x nB}
  wClusterOnDB <- Matrix::t(x = {Matrix::t(x = clusterOnDB2) * {finestAreas}})

  # group finest resolution into clusters {nC x nGibbs}
  if (is(object = yFinest, class2 = "ff_matrix")) {
    i1 <- NULL
    i2 <- NULL

    yAv <- ff::ffcolapply(EXPR = wClusterOnDB %*% yFinest[,i1:i2],
                          X = yFinest, RETURN = TRUE, CFUN = "ccbind",
                          FF_RETURN = FALSE)
  } else {
    yAv <- wClusterOnDB %*% yFinest
  }

  # {nB x nGibbs}
  # repeat rows to put back into nB dimension
  yAv <- yAv[idx,]

  # {nGibbs x nB}
  if (is(object = yAv, class2 = "Matrix")) {
    yAv <- as.matrix(x = Matrix::t(x = yAv))
  } else {
    yAv <- t(x = yAv)
  }

  # {nB x nSource}
  finestOnSource <- t(x = finestOnSource)

  # {nGibbs x nSource}
  yAv <- yAv %*% finestOnSource

  # {nSource}

  if (is(object = yFinest, class2 = "ff_matrix")) {
    i1 <- NULL
    i2 <- NULL
    tt <- as.numeric(x = ff::ffcolapply(EXPR = colMeans({ySource[,i1:i2] - yAv[,i1:i2]}^2),
                                        X = ySource, RETURN = TRUE, CFUN = "c",
                                        FF_RETURN = FALSE))
  } else {
    tt <- as.numeric(x = Matrix::colMeans(x = {ySource - yAv}^2))
  }

  # {nSource}
  idx2 <- drop(x = idx %*% finestOnSource)

  # {nC x 2}
  agg <- aggregate(x = tt, by = list(idx2), FUN = sum)

  agg1 <- numeric(length = nC)
  agg1[agg[,1L]] <- agg[,2L]

  # total area for each cluster {nCluster}
  totalSourceAreas <- drop(x = clusterOnDB %*% finestOnSource %*% sourceAreas)

  totalSourceAreas[totalSourceAreas <= 1e-8] <- 1.0

  agg1 <- agg1 / totalSourceAreas

  return( as.matrix(x = agg1) )

}
