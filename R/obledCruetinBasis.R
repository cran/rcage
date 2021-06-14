# Obled-Cruetin Basis Function
#
# function is not exported
#
# Performs a reweighting of the basis function so that they are orthonormal.
#
# @param ... Ignored. Used only to require named input.
#
# @param spatialData A SpatialPoints object, SpatialPolygons object, or a 
#   list of said objects. The source support.
#
# @param dB A SpatialPolygons object, integer, or NULL. The finest resolution.
#
# @param gbfObj A GBFObj object. The basis function information. 
#
# @param nw An integer object or NULL. The number of MC replicates to 
#   generate for estimating the O-C eigenfunctions. If NULL, 'spatialData' 
#   must be or include SpatialPoints data.
#
# @param localCluster A cluster object or NULL. A cluster if using parallel.
#
# @param verify A logical. If TRUE, verify that the scaling parameter is ok.
#
# @return A list object. The {n x r} matrix of basis functions and the
#  {r x r} OC weighting matrix
#
#' @include choleskyInvW.R generatingBasis.R
#
.obledCruetinBasis <- function(..., 
                               spatialData, 
                               dB,  
                               gbfObj,  
                               nw,  
                               localCluster,
                               verify = FALSE) {

  message("calculating OC basis functions")

  # PwLambdaW^{-1/2} { r x r }
  cholInvW <- .cholesky_W(spatialData = spatialData,
                          dB = dB,
                          gbfObj = gbfObj,
                          nw = nw,
                          verify = verify)

  # reset bandwidth in case it was modified by .cholesky_W
  gbfObj@args$w <- cholInvW$w

  # g(s) { nSpatial x r }
  basis <- .generatingBasis(spatialData = spatialData,
                            nw = nw,
                            gbfObj = gbfObj,
                            db = min(10000L, nw),
                            localCluster = localCluster)

  ## normalize the phiOC basis

  phiOC <- basis %*% cholInvW$matrix

  normal <- t(x = phiOC) %*% phiOC

  dnorm <- diag(x = normal)

  for (i in 1L:length(x = dnorm)) {

    if (abs(dnorm[i]) <= 1e-8) next

    phiOC[,i] <- phiOC[,i] / sqrt(x = dnorm[i])

  }

  return( list("basis" = basis, 
               "OCnorm" = cholInvW$matrix,
               "phiOC" = phiOC,
               "gbfObj" = gbfObj) )

}
