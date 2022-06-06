# Verify that Covariate Input is Appropriate in Content and Structure
#
# method is not exported
#
# Checks covariate is numeric and of sufficient length. Returns a matrix.
#
# @param spatialData A SpatialPolygons object, SpatialPoints object, or a
#   list of said objects. The source support.
#
# @param x A matrix, numeric, data.frame, list, or character
#   object. The covariates.
#
# @param ... ignored
#
# @return the covariate matrix
#
# @name verifyCovariates
#
setGeneric(name = ".verifyCovariates",
           def = function(spatialData, x, ...) { 
               standardGeneric(".verifyCovariates") 
             })

# default method returns error
setMethod(f = ".verifyCovariates",
          signature = c(spatialData = "ANY",
                        x = "ANY"),
          definition = function(spatialData, x, ...) { stop("not allowed") })

# ensure provided matrix is appropriate dimensions
setMethod(f = ".verifyCovariates",
          signature = c(spatialData = "Spatial",
                        x = "matrix"),
          definition = function(spatialData, x, ...) {

              # ensure the length of x matches the length of spatialData
              if (nrow(x = x) != length(x = spatialData)) {
                stop('x must be provided on dimension of spatialData',
                     call. = FALSE)
              }

              if (any(is.na(x = x))) stop('x must be complete')

              return( x )
            })


# if provided as vector convert to matrix and call matrix method
setMethod(f = ".verifyCovariates",
          signature = c(spatialData = "Spatial",
                        x = "numeric"),
          definition = function(spatialData, x, ...) {

              return( .verifyCovariates(spatialData = spatialData, 
                                        x =  matrix(data = x, ncol = 1L), ...) )

            })

# if data.frame provided convert to data.matrix and call matrix method
setMethod(f = ".verifyCovariates",
          signature = c(spatialData = "Spatial",
                        x = "data.frame"),
          definition = function(spatialData, x, ...) {

              return( .verifyCovariates(spatialData = spatialData,
                                        x = unname(obj = data.matrix(frame = x)), 
                                        ...) )

            })

# if character provided extract from @data and call matrix method
setMethod(f = ".verifyCovariates",
          signature = c(spatialData = "Spatial",
                        x = "character"),
          definition = function(spatialData, x, ...) {

              if (!all(x %in% colnames(x = spatialData@data))) {
                stop("x not in Spatial object data.frame", call. = FALSE)
              }

              x <- spatialData@data[, x, drop = FALSE]

              return( .verifyCovariates(spatialData = spatialData,
                                        x = x, ...) )
            })

# if NULL provided set covariates as intercept only
setMethod(f = ".verifyCovariates",
          signature = c(spatialData = "Spatial",
                        x = "NULL"),
          definition = function(spatialData, x, ...) {

              message("intercept only model")

              return( matrix(data = 1.0, 
                             nrow = nrow(x = spatialData@data), 
                             ncol = 1L) )

            })

# if matrix provided for multiple Spatial objects verify dimension
setMethod(f = ".verifyCovariates",
          signature = c(spatialData = "list",
                        x = "matrix"),
          definition = function(spatialData, x, ...) {

              n <- 0L

              for (i in 1L:length(x = spatialData)) {

                n <- n + length(x = spatialData[[ i ]])

              }

              if (n != nrow(x = x)) stop("x dim error", call. = FALSE)

              if (any(is.na(x = x))) stop('x must be complete')

              return( x )

            })

# if numeric provided for multiple Spatial objects convert to matrix and
#   call matrix method
setMethod(f = ".verifyCovariates",
          signature = c(spatialData = "list",
                        x = "numeric"),
          definition = function(spatialData, x, ...) {

              return( .verifyCovariates(spatialData = spatialData,
                                        x = matrix(data = x, ncol = 1L), ...) )

            })

# if character provided for multiple Spatial objects call character method
#   for each Spatial object
setMethod(f = ".verifyCovariates",
          signature = c(spatialData = "list",
                        x = "character"),
          definition = function(spatialData, x, ...) {

              result <- NULL

              for (i in 1L:length(x = spatialData)) {

                xt <- .verifyCovariates(spatialData = spatialData[[ i ]],
                                        x = x, ...)

                result <- rbind(result, xt)

              }

              return( result )

            })

# if list provided for multiple Spatial objects call appropriate method for
#  each list element
setMethod(f = ".verifyCovariates",
          signature = c(spatialData = "list",
                        x = "list"),
          definition = function(spatialData, x, ...) {

              result <- NULL

              for (i in 1L:length(x = spatialData)) {

                xt <- .verifyCovariates(spatialData = spatialData[[ i ]],
                                        x = x[[ i ]], ...)

                result <- rbind(result, xt)

              }

              return( result )

            })

# if NULL provided set covariates as intercept only
setMethod(f = ".verifyCovariates",
          signature = c(spatialData = "list",
                        x = "NULL"),
          definition = function(spatialData, x, ...) {

              message("intercept only model\n")

              n <- 0L

              for (i in 1L:length(x = spatialData)) {

                n <- n + length(x = spatialData[[ i ]])

              }

              return( matrix(data = 1.0, nrow = n, ncol = 1L) )

            })
