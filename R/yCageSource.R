setGeneric(name = ".yCageSourcePrep",
           def = function(spatialData, ...) { 
                   standardGeneric(".yCageSourcePrep") 
                 })

setMethod(f = ".yCageSourcePrep",
          signature = c(spatialData = "ANY"),
          definition = function(spatialData, ...) { stop("not allowed") })

# if CAGE and source is points only, extract points component from
# Gibbs (diagonal matrix)
# {nSource x nSource}
setMethod(f = ".yCageSourcePrep",
          signature = c(spatialData = "SpatialPoints"),
          definition = function(spatialData, ...) {

              indicator <- diag(nrow = length(x = spatialData), 
                                ncol = length(x = spatialData))

              areas <- rep(x = 1.0, times = length(x = spatialData))

              return( list("indicator" = indicator,
                           "areas" = areas) )
            })



# if CAGE and multi-resolution source support, extract points component from
# Gibbs (block diagonal matrix)
setMethod(f = ".yCageSourcePrep",
          signature = c(spatialData = "list"),
          definition = function(spatialData, ...) {

               # total number of source supports
               n <- sum(sapply(X = spatialData, FUN = length))

               n1 <- 0L
               for (i in 1L:length(x = spatialData)) {

                 n2 <- n1 + length(x = spatialData[[ i ]])

                 if (is(object = spatialData[[ i ]], class2 = "SpatialPoints")) {

                   n21 <- n2 - n1

                   indicator <- matrix(data = 0.0, 
                                       nrow = n21, 
                                       ncol = n)

                   indicator[,{n1+1L}:n2] <- diag(nrow = n21, ncol = n21)

                   areas <- rep(x = 1.0, times = n21)
                   break
                 }
 
                 n1 <- n2
               }

               return( list("indicator" = indicator,
                            "areas" = areas) )
             })

setGeneric(name = ".yDCAGESourcePrep",
           def = function(spatialData, ...) { 
                   standardGeneric(".yDCAGESourcePrep") 
                 })

setMethod(f = ".yDCAGESourcePrep",
          signature = c(spatialData = "ANY"),
          definition = function(spatialData, ...) { stop("not allowed") })


# if DCAGE and source is polygons only, extract polygons component from
# Gibbs (diagonal matrix)
# {nSource x nSource} which is equivalent to {nSource x nSource}
setMethod(f = ".yDCAGESourcePrep",
          signature = c(spatialData = "SpatialPolygons"),
          definition = function(spatialData, ...) {

              indicator <- diag(nrow = length(x = spatialData), 
                                ncol = length(x = spatialData))

              areas <- .getArea(spatialData = spatialData, byid = TRUE)

              return( list("indicator" = indicator,
                           "areas" = areas) )
            })



# if DCAGE and multi-resolution source support, extract polygons component from
# Gibbs (block diagonal matrix)
setMethod(f = ".yDCAGESourcePrep",
          signature = c(spatialData = "list"),
          definition = function(spatialData, ...) {

               # total number of source supports
               n <- sum(sapply(X = spatialData, FUN = length))

               n1 <- 0L
               for (i in 1L:length(x = spatialData)) {

                 n2 <- n1 + length(x = spatialData[[ i ]])

                 if (is(object = spatialData[[ i ]], 
                        class2 = "SpatialPolygons")) {

                   n21 <- n2 - n1

                   indicator <- matrix(data = 0.0, 
                                       nrow = n21, 
                                       ncol = n)

                   indicator[,{n1+1L}:n2] <- diag(nrow = n21, ncol = n21)

                   areas <- .getArea(spatialData = spatialData[[ i ]], byid = TRUE)
                   break
                 }
 
                 n1 <- n2
               }

               return( list("indicator" = indicator,
                            "areas" = areas) )
             })


setGeneric(name = ".yCageSource",
           def = function(ffdir, ...) { standardGeneric(".yCageSource") })

setMethod(f = ".yCageSource",
          signature = c(ffdir = "ANY"),
          definition = function(ffdir, ...) { stop("not allowed") })

setMethod(f = ".yCageSource",
          signature = c(ffdir = "NULL"),
          definition = function(ffdir, 
                                ..., 
                                gibbsSamples,
                                x,
                                spatialData,
                                cage,
                                basis,
                                H) { 

              if (cage) {
                prep <- .yCageSourcePrep(spatialData = spatialData)
              } else {
                prep <- .yDCAGESourcePrep(spatialData = spatialData)
              }

              prep[[ "ySource" ]] <- t(x = prep$indicator %*% 
                                           {x %*% gibbsSamples$beta +
                                           H %*% gibbsSamples$xi + 
                                           basis$phiOC %*% gibbsSamples$eta})

              return( prep )
            })

setMethod(f = ".yCageSource",
          signature = c(ffdir = "character"),
          definition = function(ffdir, 
                                ..., 
                                gibbsSamples,
                                x,
                                spatialData,
                                cage,
                                basis,
                                H) { 

              if (cage) {
                prep <- .yCageSourcePrep(spatialData = spatialData)
              } else {
                prep <- .yDCAGESourcePrep(spatialData = spatialData)
              }
              i1 <- NULL
              i2 <- NULL

              prep[[ "ySource" ]] <- ff::ffcolapply(EXPR = prep$indicator %*% 
                                                           {x %*% gibbsSamples$beta[,i1:i2] +
                                                            H %*% gibbsSamples$xi[,i1:i2] + 
                                                            basis$phiOC %*% gibbsSamples$eta[,i1:i2]},
                                                     X = gibbsSamples$beta, 
                                                     RETURN = TRUE, 
                                                     FF_RETURN = TRUE,
                                        RETROW = nrow(x = prep$indicator))

             prep[[ "ySource" ]] <- ff::vt(x = prep[[ "ySource" ]])

              return( prep )
            })
