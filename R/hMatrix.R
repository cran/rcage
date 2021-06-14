#  Calculate Change of Support matrix
#
#  method is not exported
#
#  @param spatialData A SpatialPoints object, SpatialPolygons object, or 
#    list of said objects. {nSpatial}
#
#  @param dB A SpatialPolygons object, integer, or NULL. The finest resolution. {nB}
#
setGeneric(name = ".hMatrix",
           def = function(spatialData, dB, ...) { standardGeneric(".hMatrix") })

# if combination of inputs is not explicitly defined the combination is
# forbidden
setMethod(f = ".hMatrix",
          signature = c(spatialData = "ANY",
                        dB = "ANY"),
          definition = function(spatialData, dB, ...) { stop("not allowed") })

# Spatial and NULL indicates that dB is the spatialData;
# change of support is diagonal
setMethod(f = ".hMatrix",
          signature = c(spatialData = "Spatial",
                        dB = "NULL"),
          definition = function(spatialData, dB, ...) {
              ns <- length(x = spatialData)
              return( diag(x = 1.0, nrow = ns, ncol = ns) )
            })

#' @importFrom rgeos gCovers
# SpatialPoints and SpatialPolygons indicates that dB is not the spatialData
# obtain set membership to determine change of support
setMethod(f = ".hMatrix",
          signature = c(spatialData = "SpatialPoints",
                        dB = "SpatialPolygons"),
          definition = function(spatialData, dB, ...) {

              # {nSpatial x nB}
              h <- rgeos::gCovers(spgeom1 = dB, 
                                  spgeom2 = spatialData,  
                                  byid = TRUE)

              # point data lies in more than one grid
              # this happens if a point is on the boundary
              # randomly assign this point to be in only one grid
              moreThanOne <- rowSums(x = h) > 1.0
              if (any(moreThanOne)) {
                for (i in 1L:length(x = moreThanOne)) {
                  if (!moreThanOne[i]) next
                  it <- sample(x = which(h[i,] > 0), size = 1L)
                  h[i,] <- 0.0
                  h[i,it] <- 1.0
                }
              }

              # ensure that matrix returned in dimensions expected
              # {nB x nSpatial}
              if (nrow(x = h) != length(x = dB)) {
                if( {ncol(x = h) == length(x = dB)} &&
                    {nrow(x = h) == length(x = spatialData)} ) {
                  h <- t(x = h)
                } else {
                  stop("dim error in h", call. = FALSE)
                }
              }

              # {nB x nPts}
              return( h )

            })

#' @importFrom rgeos gIntersects gTouches gArea
#' @importFrom raster intersect
# SpatialPolygons and SpatialPolygons indicates that dB is not the spatialData
# obtain overlap to determine change of support
setMethod(f = ".hMatrix",
          signature = c(spatialData = "SpatialPolygons",
                        dB = "SpatialPolygons"),
          definition = function(spatialData, dB, ...) {

              # {nPoly x nB}
              h <- rgeos::gIntersects(spgeom1 = dB, 
                                      spgeom2 = spatialData,  
                                      byid = TRUE) -
                   rgeos::gTouches(spgeom1 = dB, 
                                   spgeom2 = spatialData,  
                                   byid = TRUE)

              den <- rgeos::gArea(spgeom = sp::SpatialPolygons(Srl = spatialData@polygons), 
                                  byid = TRUE)

              for (i in 1L:nrow(x = h)) {
                it <- which(h[i,] > 0)
                if (length(x = it) == 0L) next
                for (j in it) {
                  ras <- raster::intersect(x = dB[j,], y = spatialData[i,])
                  num <- rgeos::gArea(spgeom = sp::SpatialPolygons(Srl = ras@polygons))

                  h[i,j] <- num/den[i]

                }
              }

              # if partitioning is complete and correctly defined,
              # the total should be exactly 1 for any polygon covered by dB
              tt <- rowSums(x = h)
              if (any(tt > 0.0 & tt < 1.0)) {
                warning("finest resolution grid does not appear to partition the source")
              }

              # {nB x nPoly}
              return( t(x = h) )
            })


# multi-resolution data and finest resolution is one of the source supports
setMethod(f = ".hMatrix",
          signature = c(spatialData = "list",
                        dB = "integer"),
          definition = function(spatialData, dB, ...) {

              h <- matrix(data = 0.0, 
                          nrow = length(x = spatialData[[ dB ]]), 
                          ncol = 0L)

              for (i in 1L:length(x = spatialData)) {
                if (i == dB) {
                  # when itself, diagonal
                  h1 <- .hMatrix(spatialData = spatialData[[ i ]],
                                 dB = NULL)
                } else {
                  # for all others, find overlap
                  h1 <- .hMatrix(spatialData = spatialData[[ i ]],
                                 dB = spatialData[[ dB ]])
                }
                # {nB x nSpatial}
                h <- cbind(h, h1)
              }

              return( h )

            })

# multi-resolution data and finest resolution is polygon not in source
setMethod(f = ".hMatrix",
          signature = c(spatialData = "list",
                        dB = "SpatialPolygons"),
          definition = function(spatialData, dB, ...) {

              h <- matrix(data = 0.0, nrow = length(x = dB), ncol = 0L)

              for (i in 1L:length(x = spatialData)) {

                h1 <- .hMatrix(spatialData = spatialData[[ i ]],
                               dB = dB)

                # {nB x nSpatial}
                h <- cbind(h, h1)
              }

              return( h )

            })
