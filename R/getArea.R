#  Calculate Area of Spatial object(s)
#
#  method is not exported
#
#  @param spatialData A SpatialPoints or SpatialPolygons object.
#
#  @param ... ignored
#
#  @param byid A logical object. TRUE indicates that areas for each
#    areal unit are returned. FALSE indicates total area.
#
#  @return areas
#
#' @import methods
#' @import sp
setGeneric(name = ".getArea", 
           def = function(spatialData, ...) { standardGeneric(".getArea") })

# default method returns an error
setMethod(f = ".getArea",
          signature = c(spatialData = "ANY"),
          definition = function(spatialData, ...) { stop("not allowed") })

# if SpatialPolygons, use rgeos methods to calculate area of each polygon
#
#' @importFrom rgeos gArea
#' @import sp
setMethod(f = ".getArea",
          signature = c(spatialData = "SpatialPolygons"),
          definition = function(spatialData, ..., byid) {

              spObj <- sp::SpatialPolygons(Srl = spatialData@polygons)

              areas <- tryCatch(expr = rgeos::gArea(spgeom = spObj, 
                                                    byid = byid),
                                error = function(e){
                                          stop("unable to determine areas\n",
                                               e$message, call. = FALSE)
                                        })

              return( unname(obj = areas) )

            })

# if SpatialPoints, use rgeos methods to calculate hull and call 
#   SpatialPolygons method
#
#' @importFrom rgeos gConvexHull
setMethod(f = ".getArea",
          signature = c(spatialData = "SpatialPoints"),
          definition = function(spatialData, ..., byid) { 

              if (byid) return( rep(x = 1.0, times = length(x = spatialData)) )

              hull <- tryCatch(expr = rgeos::gConvexHull(spgeom = spatialData),
                                error = function(e){
                                          stop("unable to determine convex hull\n",
                                               e$message, call. = FALSE)
                                        })

              return( .getArea(spatialData = hull, byid = FALSE) )

            })
