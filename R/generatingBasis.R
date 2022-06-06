#  Generating Basis Functions
#
#  method is not exported
#
#  @param spatialData A Spatial object as defined by package sp
#    currently limited to SpatialPoints and SpatialPolygons or list thereof
#
#  @param ... ignored.
#
#  @param gbfObj A GBFObj object. Information about basis functions.
#
#  @param nw A numeric object or NULL. If numeric, the number of MC replicates.
#
#  @param db A numeric object. The minimum number of samples to take at 
#    one call to sp::spsample. This is only to improve performance.
#
#  @param localCluster A cluster object or NULL. If parallel methods 
#    requested, a previously established cluster.
#
#  @return matrix of basis functions { nSpatial x r }
#
#' @include GBFObj.R
setGeneric(name = ".generatingBasis",
           def = function(spatialData, nw, ...) { 
                   standardGeneric(".generatingBasis") 
                 })

# default method returns error
setMethod(f = ".generatingBasis",
          signature = c(spatialData = "ANY",
                        nw = "ANY"),
          definition = function(spatialData, nw, ...) { stop("not allowed") })

# calculate basis on coordinates of SpatialPoints
#' @importFrom sp coordinates
setMethod(f = ".generatingBasis",
          signature = c(spatialData = "SpatialPoints",
                        nw = "ANY"),
          definition = function(spatialData, nw, ..., gbfObj) {

              message("\tcalculating basis using SpatialPoints data")

              # calculate basis function at coordinates of SpatialPoints
              return( .basis(crd = sp::coordinates(obj = spatialData),
                             gbfObj = gbfObj) )

            })

#' @importFrom sp spsample coordinates
.iteration <- function(i, ..., nw, spatialData, db, gbfObj) {

  # generate nw samples from areal unit i

  hold <- matrix(data = 0.0, nrow = 0L, ncol = 2L)

  failed <- FALSE

  msg <- paste("unable to generate sample points in area", i, "\n")

  while (nrow(x = hold) < nw) {

    # sample areal unit to obtain db samples
    points <- tryCatch(expr = sp::spsample(x = spatialData[i,],
                                           n = db,
                                           type = "random"),
                       error = function(e) {
                                 if (failed) {
                                   # if second failure, stop with error
                                   stop(msg, e$message, call. = FALSE)
                                 } else {
                                   failed <- TRUE
                                 }
                               })

    # if sampling failed, try again
    if (failed) next

    # add coordinates of sampled points to matrix; will exit when 
    # number of samples is at least nw
    hold <- rbind(hold, sp::coordinates(obj = points))

  }

  # calculates psi(s_i) i = 1, m and average over MC samples
  return( colMeans(x = .basis(crd = hold[1L:nw,,drop = FALSE],
                              gbfObj = gbfObj)) )
}

#' @importFrom parallel clusterExport parSapply
# SpatialPolygons provided use Monte Carlo
setMethod(f = ".generatingBasis",
          signature = c(spatialData = "SpatialPolygons",
                        nw = "numeric"),
          definition = function(spatialData,
                                nw, 
                                ..., 
                                gbfObj,
                                db,
                                localCluster) {

              message("\tapproximating GBF using MC sampling")

              # number of areal units
              nA <- length(x = spatialData)

              # number of basis function
              r <- nrow(x = gbfObj@args$knots)

              ## generate basis functions

              db <- min(db, nw)

              if (!is.null(x = localCluster)) {

                # if using parallel methods and user provides a basis
                # function, export to cluster

                if (!{gbfObj@name %in% c(".bisquare", ".wendland", ".gauss")}) {
                  parallel::clusterExport(cl = localCluster, 
                                          varlist = list(gbfObj@name))
                }

                # for each areal unit, estimate basis function
                basis <- t(x = parallel::parSapply(cl = localCluster, 
                                                   X = 1L:nA, 
                                                   FUN = .iteration,
                                                   nw = nw, 
                                                   spatialData = spatialData, 
                                                   db = db, 
                                                   gbfObj = gbfObj))

              } else {

                basis <- matrix(data = 0.0, nrow = nA, ncol = r)

                for (a in 1L:nA) {

                  basis[a,] <- .iteration(i = a, 
                                          nw = nw, 
                                          spatialData = spatialData, 
                                          db = db, 
                                          gbfObj = gbfObj)

                }

              }

              if ({nrow(x = basis) != nA} || {ncol(x = basis) != r}) {
                if ({ncol(x = basis) == nA} && {nrow(x = basis) == r}) {
                  basis <- t(x = basis)
                } else {
                  stop("dim error in basis", call. = FALSE)
                }
              }

              return( basis )
            })

# if a list of spatial data is provided and points data not included or are
# insufficient to estimate psi, use Monte Carlo for each component of list
setMethod(f = ".generatingBasis",
          signature = c(spatialData = "list",
                        nw = "numeric"),
          definition = function(spatialData, nw, ...) {

              basis <- NULL

              for (i in 1L:length(x = spatialData)) {

                basis <- rbind(basis, 
                               .generatingBasis(spatialData = spatialData[[ i ]], 
                                                nw = nw, ...))
              }

              return( basis )

            })

# if a list of spatial data is provided and points are sufficient to estimate
# psi for all elements, use points data for each component of list
setMethod(f = ".generatingBasis",
          signature = c(spatialData = "list",
                        nw = "NULL"),
          definition = function(spatialData, nw, ..., gbfObj) {

              # identify which element of list is the points data
              for (i in 1L:length(x = spatialData)) {
                if (is(object = spatialData[[ i ]], class2 = "SpatialPoints")) {
                  ip <- i
                  break
                }
              }

              # number of knots
              r <- nrow(x = gbfObj@args$knots)

              # coordinates of points
              crds <- sp::coordinates(obj = spatialData[[ ip ]])

              basis <- NULL

              for (i in 1L:length(x = spatialData)) {

                if (i == ip) {
                  # use points to calculate basis
                  basis <- rbind(basis, 
                                 .generatingBasis(spatialData = spatialData[[ i ]], 
                                                  nw = NULL, ...))
                } else {

                  # determine which points lie in each polygon
                  # {nPoly x nPts}
                  h <- rgeos::gContains(spgeom1 = spatialData[[ i ]], 
                                        spgeom2 = spatialData[[ ip ]],  
                                        byid = TRUE)

                  for (j in 1L:nrow(x = h)) {
                    # use subset of points that lie in the jth polygon 
                    # to calculate basis
                    tt <- colMeans(x = .basis(crd = crds[h[j,],,drop=FALSE],
                                              gbfObj = gbfObj))

                    if (length(x = tt) != r) stop("dim issue")

                    basis <- rbind(basis, tt)

                  }
                }
              }

              return( basis )

            })
