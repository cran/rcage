#  Method to Verify and/or Calculate Knots
#
#  method is not exported
#
#  @param knots A numeric or numeric matrix. If numeric, the number of 
#    knots to be generated. If matrix, the x,y coordinates of knots.
#
#  @param \dots ignored.
#
#  @param spatialData A SpatialPoints, SpatialPolygons, or list of said objects.
#    The source support.
#
#  @param dB A SpatialPolygons object or NULL. The finest resolution.
#
#  @param longlat A logical object. If TRUE, spatial data is in long/lat
#    coordinate system.
#
#  @return \{nKnots x 2\} matrix of knots ordered according to y coordinate
#
setGeneric(name = ".knots",
           def = function(knots, ...) { standardGeneric(".knots") })

# default method results in an error
setMethod(f = ".knots",
          signature = c(knots = "ANY"),
          definition = function(knots, ...) { stop("not allowed") })

# if matrix provided by user, ensure that there are two columns and
#   order the rows according to the second coordinate.
setMethod(f = ".knots",
          signature = c(knots = "matrix"),
          definition = function(knots, ...) {

              if (ncol(x = knots) != 2L) stop("knots dim error", call. = FALSE)

              return( knots[order(knots[,2L]),,drop=FALSE] )

            })

# if provided by user as a data.frame, convert to matrix and call matrix method
setMethod(f = ".knots",
          signature = c(knots = "data.frame"),
          definition = function(knots, ...) {

              knots <- unname(obj = data.matrix(frame = knots))

              return( .knots(knots = knots, ...) )

            })

# if not provided as a matrix, assume value is the number of knots to generate
#   from provided Spatial object
#
#' @importFrom fields cover.design
#' @import sp
setMethod(f = ".knots",
          signature = c(knots = "numeric"),
          definition = function(knots, ..., spatialData, dB, longlat) {

              if (length(x = knots) > 1L) {
                stop('must provide matrix of knots or # of knots to generate', 
                     call. = FALSE)
              }

              # generate samples from spatialData/dB

              if (is.null(x = dB)) {
                # if dB is NULL, using source support data as finest resolution 
                if (is(object = spatialData, class2 = "SpatialPoints")) {
                  # get convex hull around points as sampling region
                  message("sampling convex hull of spatialData to obtain knots")
                  dB <- rgeos::gConvexHull(spgeom = spatialData)
                } else {
                  message("sampling spatialData to obtain knots")
                  dB <- spatialData
                }
              } else if (is.numeric(x = dB)) {
                # if dB is integer, using a component of source support data 
                # as finest resolution.
                if (is(object = spatialData[[ dB ]], class2 = "SpatialPoints")) {
                  message("sampling convex hull of spatialData element ", dB, 
                          " to obtain knots")
                  dB <- rgeos::gConvexHull(spgeom = spatialData[[ dB ]])
                } else {
                  message("sampling spatialData element ", dB, 
                          " to obtain knots")
                  dB <- spatialData[[ dB ]]
                }
              } else {
                message("sampling dB to obtain knots")
              }

              # sample region uniformly
              sm <- tryCatch(expr = sp::spsample(x = dB,
                                                 n = 100L*knots,
                                                 type = "random"),
                             error = function(e) {
                                       stop('unable to sample data\n',
                                            e$message, call. = FALSE)
                                     })

              # define function to calculate pairwise distance between points
              disfunc <- function(x1, x2) {
                           return(sp::spDists(x = x1, 
                                              y = x2, 
                                              longlat = longlat))
                         }

              message("using fields::cover.design() to obtain knots")

              # find the set of points on a discrete grid that minimizes a
              # geometric space-filling criterion.
              knots <- tryCatch(expr = fields::cover.design(R = sp::coordinates(obj = sm), 
                                                            nd = knots,
                                                            DIST = disfunc),
                                error = function(e) {
                                          stop('unable to obtain knots\n',
                                               e$message, call. = FALSE)
                                        })

              return( .knots(knots = knots$design, ...) )

            })
