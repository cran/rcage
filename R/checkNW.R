# Ensure that the provided nw value (used to sample spatial region when
# estimating W and Psi(A)) is appropriate.
#
# spatialData	nw	W/Psi(A) calculated using
#--------------------------------------------------------
# pts		NULL	pts data
# pts		nMC	nMC samples of convexHull of pts
# poly		NULL	20000 samples of poly
# poly		nMC	nMC samples of poly
# list		NULL	n = 1; pts element
#			n = 0; 20000 samples of dB
# list		nMC	nMC samples of dB
setGeneric(name = ".checkNW",
           def = function(spatialData, nw, ...) { standardGeneric(".checkNW") })

# if not expressly allowed, input combination is forbidden
setMethod(f = ".checkNW",
          signature = c(spatialData = "ANY",
                        nw = "ANY"),
          definition = function(spatialData, nw, ...) { stop("not allowed") })

# SpatialPoints is provided and is sufficient to estimate W and Psi
setMethod(f = ".checkNW",
          signature = c(spatialData = "SpatialPoints",
                        nw = "NULL"),
          definition = function(spatialData, nw, ...) {

              message("D_S is assumed to provide sufficient coverage ",
                      "for direct calculation of W and psi(A)")

              return( nw )
            })

# SpatialPoints is provided but is not sufficient to estimate W and 
# psi(A)
setMethod(f = ".checkNW",
          signature = c(spatialData = "SpatialPoints",
                        nw = "numeric"),
          definition = function(spatialData, nw, ...) {

              nw <- as.integer(x = round(x = nw, digits = 0L))

              if (nw == 0L) {

                nw <- .checkNW(spatialData = spatialData, nw = NULL)

              } else {

                if (nw < 1000L) {
                  message("nw must be >= 1000 when D_S ",
                          "does not provide sufficient coverage for the ",
                          "direct calculation W and psi(A); ",
                          "nw set to 20000")
                  nw <- 20000L
                }

                message("W and psi(A) will be estimated using ", 
                        nw, " MC samples")
              }

              return( nw )

            })

# SpatialPolygons provided without SpatialPoints with sufficient coverage but
# nw was not set by user
setMethod(f = ".checkNW",
          signature = c(spatialData = "SpatialPolygons",
                        nw = "NULL"),
          definition = function(spatialData, nw, ...) {

              message("nw must be >= 1000 when D_S is not provided; ",
                      "nw set to 20000")
              nw <- 20000L

              message("W and psi(A) will be estimated using ", 
                      nw, " MC samples")

              return( nw )

            })

# SpatialPolygons provided without SpatialPoints with sufficient coverage
# nw was set by user
setMethod(f = ".checkNW",
          signature = c(spatialData = "SpatialPolygons",
                        nw = "numeric"),
          definition = function(spatialData, nw, ...) {

              nw <- as.integer(x = round(x = nw, digits = 0L))

              if (nw < 1000L) {
                message("nw must be >= 1000 when D_S is not provided; ",
                        "nw set to 20000")
                nw <- 20000L
              }

              message("W and psi(A) will be estimated using ", 
                      nw, " MC samples")

              return( nw )

            })

# multiresolution data provided and points data is sufficient for calculating
# W and psi(A) (User specified nw as NULL)
#' @importFrom methods is
setMethod(f = ".checkNW",
          signature = c(spatialData = "list",
                        nw = "NULL"),
          definition = function(spatialData, nw, ...) {

              for (i in 1L:length(x = spatialData)) {

                if (is(object = spatialData[[ i ]], 
                       class2 = "SpatialPoints")) {

                  return( .checkNW(spatialData = spatialData[[ i ]], 
                                   nw = NULL) )

                }
              }

              return( .checkNW(spatialData = spatialData[[ 1L ]], nw = NULL) )

            })

# multiresolution data provided and points data is insufficient for calculating
# W and Psi. (User specified nw)
setMethod(f = ".checkNW",
          signature = c(spatialData = "list",
                        nw = "numeric"),
          definition = function(spatialData, nw, ...) {

              nw <- as.integer(x = round(x = nw, digits = 0L))

              if (nw == 0L) {

                return( .checkNW(spatialData = spatialData, nw = NULL) )

              } else {

                if (nw < 1000L) {
                  message("nw must be >= 1000 when D_S is not provided or ",
                          "does not provide sufficient coverage for the ",
                          "direct calculation W and psi(A); ",
                          "nw set to 20000")
                  nw <- 20000L
                }

                message("W and psi(A) will be estimated using ", 
                        nw, " MC samples")
              }

              return( nw )

            })
