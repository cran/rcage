# Ensure that spatialData takes expected form. Currently limited to
# SpatialPoints, SpatialPointsDataFrame, SpatialPolygons, 
# SpatialPolygonsDataFrame or list of said objects
# 
#' @import sp
setGeneric(name = ".checkSpatialData",
           def = function(spatialData, ...) { 
                   standardGeneric(".checkSpatialData") 
                 })

# If not expressly defined, input is forbidden
setMethod(f = ".checkSpatialData",
          signature = c(spatialData = "ANY"),
          definition = function(spatialData, ...) { 
                         stop("bad spatialData", call. = FALSE) 
                       })

# SpatialPointsDataFrame provided -- ok
setMethod(f = ".checkSpatialData",
          signature = c(spatialData = "SpatialPointsDataFrame"),
          definition = function(spatialData, ...) { return( spatialData )})

# SpatialPolygonsDataFrame provided -- ok
setMethod(f = ".checkSpatialData",
          signature = c(spatialData = "SpatialPolygonsDataFrame"),
          definition = function(spatialData, ...) { return( spatialData )})

# Convert SpatialPoints to SpatialPointsDataFrame
setMethod(f = ".checkSpatialData",
          signature = c(spatialData = "SpatialPoints"),
          definition = function(spatialData, ...) { 

             df <- data.frame("X" = rep(x = 1.0, 
                                        times = length(x = spatialData)))

             return( sp::SpatialPointsDataFrame(coords = spatialData,
                                                data = df) )
            })

# Convert SpatialPolygons to SpatialPolygonsDataFrame
setMethod(f = ".checkSpatialData",
          signature = c(spatialData = "SpatialPolygons"),
          definition = function(spatialData, ...) {

              df <- data.frame("X" = rep(x = 1.0, 
                                         times = length(x = spatialData)))

              return( sp::SpatialPolygonsDataFrame(Sr = spatialData, 
                                                   data = df) )
            })

#' @importFrom methods is
setMethod(f = ".checkSpatialData",
          signature = c(spatialData = "list"),
          definition = function(spatialData, ...) {

              for (i in 1L:length(x = spatialData)) {

                if (!is(object = spatialData[[ i ]], 
                        class2 = "SpatialPoints") &&
                    !is(object = spatialData[[ i ]], 
                        class2 = "SpatialPolygons")) {
                  stop("unrecognized object in spatialData", call. = FALSE)
                }

                spatialData[[ i ]] <- .checkSpatialData(spatialData = spatialData[[ i ]])

              }
              return( spatialData )
            })
