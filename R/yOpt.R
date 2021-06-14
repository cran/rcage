# Area Weighted Value of Interest in Each Cluster for All Gibbs Samples
#
# method is not exported
#
# @param y A matrix or ff_matrix object. The modeled outcome on dB for 
#   all kept Gibbs. {nB x nKept}
#
# @param idxit A vector object. The cluster ids for each areal unit in dB.
#    {nB}
#
# @param areas A vector object. The area of each areal unit in dB. {nB}
#
# @return A matrix or ff_matrix object. {nC x nKept}
#
setGeneric(name = ".yOpt",
           def = function(y, ...) { standardGeneric(".yOpt") })

# default method results in error
setMethod(f = ".yOpt",
          signature = c(y = "ANY"),
          definition = function(y, ...) { stop("not allowed") })

#' @importFrom stats aggregate
setMethod(f = ".yOpt",
          signature = c(y = "matrix"),
          definition = function(y, ..., idxit, areas) {

              # identify unique cluster ids {nClusters}
              unq <- sort(x = unique(x = idxit))

              # number of clusters
              nClusters <- length(x = unq)

              message(nClusters, " clusters")

              # cluster id of each areal unit
              idx <- match(x = idxit, table = unq)

              # total area for each cluster {nCluster}
              totalAreas <- stats::aggregate(x = areas,
                                             by = list(idx),
                                             FUN = sum)[,-1L]

              # area weighted value of interest for each cluster and sample
              # {nCluster x nGibbs}

              YOpt <- stats::aggregate(x = y*areas,
                                       by = list(idx),
                                       FUN = sum)[,-1L] / totalAreas
              YOpt <- unname(obj = data.matrix(frame = YOpt))

              return( YOpt )
            })

#' @import ff
#' @importFrom stats aggregate
setMethod(f = ".yOpt",
          signature = c(y = "ff_matrix"),
          definition = function(y, ..., idxit, areas) {

              # identify unique cluster ids {nClusters}
              unq <- sort(x = unique(x = idxit))

              # number of clusters
              nClusters <- length(x = unq)

              message(nClusters, " clusters")

              # cluster id of each areal unit
              idx <- match(x = idxit, table = unq)

              # total area for each cluster {nCluster}
              totalAreas <- stats::aggregate(x = areas,
                                             by = list(idx),
                                             FUN = sum)[,-1L]

              # area weighted value of interest for each cluster and sample
              # {nCluster x nGibbs}
              i1 <- NULL
              i2 <- NULL
              YOpt <- ff::ffcolapply(EXPR = data.matrix(frame = stats::aggregate(x = y[,i1:i2]*areas,
                                                                                 by = list(idx),
                                                                                 FUN = sum)[,-1L] / totalAreas),
                                     X = y, 
                                     RETURN = TRUE, 
                                     FF_RETURN = TRUE, 
                                     RETROW = nClusters)

              return( YOpt )
            })
