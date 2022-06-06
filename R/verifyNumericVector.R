# Verify that a Vector Input is Appropriate in Content and Structure
#
# method is not exported
#
# Checks that a vector is numeric and is of sufficient length. Functions
#   return the vector
#
# @param spatialData A SpatialPolygons object, SpatialPoints object,
#   or a list of said objects. The source support.
#
# @param numVec A matrix, numeric, data.frame, list, or 
#   character object. The variable of interest. 
#
# @param ... ignored.
#
# @return the vector
#
# @name verifyNumericVector
#
setGeneric(name = ".verifyNumericVector",
           def = function(spatialData, numVec, ...) { 
               standardGeneric(".verifyNumericVector") 
             })

# default method results in an error
setMethod(f = ".verifyNumericVector",
          signature = c(spatialData = "ANY",
                        numVec = "ANY"),
          definition = function(spatialData, numVec, ...) { 
              stop("not allowed") 
            })

# If a single Spatial object and matrix provided, ensure dimensions are
# appropriate
setMethod(f = ".verifyNumericVector",
          signature = c(spatialData = "Spatial",
                        numVec = "matrix"),
          definition = function(spatialData, numVec, ..., nm) {

              if (ncol(x = numVec) == 1L) {
                return( .verifyNumericVector(spatialData = spatialData,
                                             numVec = drop(numVec), ...,
                                             nm = nm) )
              }

              # check if numVec is numeric (this includes integer)
              if (!is.numeric(x = numVec[,1L])) {
                stop(nm, " is not numeric", call. = FALSE)
              }

              # ensure the length of numVec matches the length of spatialData
              if (nrow(x = numVec) != length(x = spatialData) ||
                  ncol(x = numVec) != length(x = spatialData)) {
                stop(nm, ' must be provided on dimension of spatialData',
                     call. = FALSE)
              }

              return( numVec )

            })

setMethod(f = ".verifyNumericVector",
          signature = c(spatialData = "Spatial",
                        numVec = "numeric"),
          definition = function(spatialData, numVec, ..., nm) {

              # ensure the length of numVec matches the length of spatialData
              if (length(x = numVec) != length(x = spatialData)) {
                stop(nm, ' must be provided on dimension of spatialData',
                     call. = FALSE)
              }

              return( numVec )


            })

setMethod(f = ".verifyNumericVector",
          signature = c(spatialData = "Spatial",
                        numVec = "character"),
          definition = function(spatialData, numVec, ..., nm) {

              if (length(x = numVec) != 1L) {
                stop("too many variables specified", call. = FALSE)
              }

              if (!{numVec %in% colnames(x = spatialData@data)}) {
                stop(nm, " not in spatial object's data.frame",
                     call. = FALSE)
              }

              numVec <- data.matrix(frame = spatialData@data[, numVec])

              return( .verifyNumericVector(spatialData = spatialData, 
                                           numVec = numVec,
                                           nm = nm) )

            })

setMethod(f = ".verifyNumericVector",
          signature = c(spatialData = "Spatial",
                        numVec = "data.frame"),
          definition = function(spatialData, numVec, ...) {

              return( .verifyNumericVector(spatialData = spatialData,
                                           numVec = data.matrix(frame = numVec), 
                                           ...) )

            })

setMethod(f = ".verifyNumericVector",
          signature = c(spatialData = "list",
                        numVec = "matrix"),
          definition = function(spatialData, numVec, ..., nm) {


              # count total number of data points/areas
              n <- 0L

              for (i in 1L:length(x = spatialData)) {
                n <- n + length(x = spatialData[[ i ]])
              }

              if (n != nrow(x = numVec) || n != ncol(x = numVec)) {
                stop(nm, " dim error", call. = FALSE)
              }

              return( numVec )

            })

setMethod(f = ".verifyNumericVector",
          signature = c(spatialData = "list",
                        numVec = "numeric"),
          definition = function(spatialData, numVec, ...) {

              return( .verifyNumericVector(spatialData = spatialData,
                                           numVec = matrix(data = numVec, 
                                                           ncol = 1L), ...) )

            })

setMethod(f = ".verifyNumericVector",
          signature = c(spatialData = "list",
                        numVec = "data.frame"),
          definition = function(spatialData, numVec, ...) {

              return( .verifyNumericVector(spatialData = spatialData,
                                           numVec = data.matrix(frame = numVec), 
                                           ...) )

            })

setMethod(f = ".verifyNumericVector",
          signature = c(spatialData = "list",
                        numVec = "list"),
          definition = function(spatialData, numVec, ..., nm) {

              if (length(x = numVec) != length(x = spatialData)) {
                stop(nm, " dim error", call. = FALSE)
              }

              result <- NULL

              for (i in 1L:length(x = spatialData)) {

                resp <- .verifyNumericVector(spatialData = spatialData[[ i ]],
                                             numVec = numVec[[ i ]], 
                                             nm = nm, ...)

                result <- c(result, resp)

              }

              return( result )

            })

setMethod(f = ".verifyNumericVector",
          signature = c(spatialData = "list",
                        numVec = "character"),
          definition = function(spatialData, numVec, ..., nm) {

              result <- NULL

              for (i in 1L:length(x = spatialData)) {

                resp <- .verifyNumericVector(spatialData = spatialData[[ i ]],
                                             numVec = numVec,
                                             nm = nm, ...)

                result <- c(result, resp)

              }

              return( result )

            })
