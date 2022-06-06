# Class containing all necessary basis information
#
# Class is not exported.
#
# @slot name the function name to be called to calculate basis
# @slot args a list of the non-coordinate arguments (knots, w, longlat)
#
setClass(Class = "GBFObj",
         slots = c("name" = "character",
                   "args" = "list"))

# Confirm Input gbf
#
# method is not exported
#
# Verify that user specified function for radial basis functions contains
#   correct input arguments
#
# @param gbf A "function" or "character" object. The function or function name 
#   to be used for calculating the radial basis functions
#
# @param ... ignored.
#
# @param weight A numeric object. The scaling factor or bandwidth.
#
# @param knots An integer or matrix object. The number of knots to generate
#   or a matrix of the coordinates of the knots.
#
# @param spatialData A SpatialPoints object, SpatialPolygons object, or 
#   list of said objects. The source support.
#
# @param longlat A logical object. TRUE indicates spatialData is in long/lat
#    coordinates
#
# @param dB A SpatialPolygons object, integer, or NULL. The finest 
#    resolution data.
#
# @return the GBFObj object
#
setGeneric(name = ".verifyGBF",
           def = function(gbf, ...) { standardGeneric(".verifyGBF") })

# default methods results in an error
setMethod(f = ".verifyGBF",
          signature = c(gbf = "ANY"),
          definition = function(gbf, ...) {
              stop("if provided, gbf must be a function or a function name",
                   call. = FALSE)
            })

# if function specified, verify structure of input arguments of function
#' @include knots.R
#' @importFrom sp spDists
setMethod(f = ".verifyGBF",
          signature = c(gbf = "character"),
          definition = function(gbf,  
                                ..., 
                                weight, 
                                knots, 
                                spatialData, 
                                longlat,
                                dB) {

              message(gbf, " is the generating basis function")

              if (gbf == "bisquare") gbf <- ".bisquare"
              if (gbf == "wendland") gbf <- ".wendland"
              if (gbf == "gaussian") gbf <- ".gauss"

              # function must use crd, knots, w, and ellipsis
              requiredfm <- c("crd", "knots", "w", "...")

              fm <- names(x = formals(fun = gbf))

              if (!all(requiredfm %in% fm)) {
                stop("gbf function must use formals crd, knots, w, and ellipsis",
                     call. = FALSE)
              }

              # verify and/or generate knots
              knots <- .knots(knots = knots, 
                              spatialData = spatialData,  
                              dB = dB,  
                              longlat = longlat)

              # if not provided, initiate weight to be 1.5*minimum distance
              # between knots
              if (is.null(x = weight)) {
                dis <- sp::spDists(x = knots, y = knots, longlat = longlat)

                weight <- 1.5 * min(dis[upper.tri(x = dis, diag = FALSE)])
                weight <- round(x = weight, digits = 4L)
                if (weight < 0.0001) weight <- 1e-3
                message("basis function scale initially set to ", weight)
              } else {
                message("basis function scale provided as ", weight)
              }


              # only a single weight value can be specified for default methods
              if (length(x = weight) > 1L && 
                  {gbf %in% c(".bisquare", ".wendland", ".gauss")}) {
                warning("more than 1 GBF scaling factor specified; using first")
                weight <- weight[1L]
              }

              return( new(Class = "GBFObj",
                          "name" = gbf,
                          "args" = list("w" = weight, 
                                        "knots" = knots,  
                                        "longlat" = longlat)) )

            })

# if user provided the entire function, strip out function name and call
#   character method
setMethod(f = ".verifyGBF",
          signature = c(gbf = "function"),
          definition = function(gbf, ...) {

              # cast function name as character and verify
              fname <- as.character(x = substitute(expr = gbf))

              return( .verifyGBF(gbf = fname, ...) )

            })

# generate radial basis functions
#
# function is not exported
#
# @param ... ignored.
#
# @param crd is a matrix of coordinates
#
# @param gbfObj an GBFObj object
#
.basis <- function(..., crd, gbfObj) {

  # add coordinates to argument list
  argList <- gbfObj@args
  argList[[ "crd" ]] <- crd

  tt <- do.call(what = gbfObj@name, args = argList) 

  return( tt )

}
