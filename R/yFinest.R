setGeneric(name = ".yFinestPrep",
           def = function(spatialData, dB, ...) { 
                    standardGeneric(".yFinestPrep") 
                 })

setMethod(f = ".yFinestPrep",
          signature = c(spatialData = "ANY",
                        dB = "ANY"),
          definition = function(spatialData, dB, ...) { stop("not allowed") })


# if source support is finest support, do not need to average
setMethod(f = ".yFinestPrep",
          signature = c(spatialData = "Spatial",
                        dB = NULL),
          definition = function(spatialData, dB, ..., fos) {

              # { nB x nSource }
              return( diag(nrow = nrow(x = fos), ncol = nrow(x = fos)) )

            })

# if dB is an integer, using an element of the source support as dB
# create indicator matrix to pull out appropriate components
setMethod(f = ".yFinestPrep",
          signature = c(spatialData = "list",
                        dB = "numeric"),
          definition = function(spatialData, dB, ..., fos) {

              n1 <- 0L
              for (i in 1L:length(x = spatialData)) {

                if (i == dB) {

                  n2 <- n1 + length(x = spatialData[[ i ]])

                  # { nB x nSource }
                  # do not need to average over source support to obtain x
                  # on dB
                  avgMatrix <- matrix(data = 0.0, 
                                      nrow = {n2 - n1}, 
                                      ncol = ncol(x = fos))

                  avgMatrix[,{n1+1L}:n2] <- diag(nrow = {n2 - n1}, 
                                                 ncol = {n2 - n1})

                  break
                }

                n1 <- n1 + length(x = spatialData[[ i ]])
              }

              return( avgMatrix )
            })

# if dB is SpatialPolygon finestOnSource projects between spaces
setMethod(f = ".yFinestPrep",
          signature = c(spatialData = "ANY",
                        dB = "SpatialPolygons"),
          definition = function(spatialData, dB, ..., fos) {

              # need to average source support {nB x nSource}
              # identify non-zero elements of finestOnSource           
              avgMatrix <- {fos > 1e-8}*1.0

              # count the number of non-zero elements for each nB
              sums <- rowSums(x = avgMatrix)
              sums[sums<1e-8] <- 1.0

              # scale by total count
              avgMatrix <- avgMatrix / sums

              return( avgMatrix )
            })

setGeneric(name = ".yFinest",
           def = function(ffdir, ...) { standardGeneric(".yFinest") })

setMethod(f = ".yFinest",
          signature = c(ffdir = "ANY"),
          definition = function(ffdir, ...) { stop("not allowed") })

setMethod(f = ".yFinest",
          signature = c(ffdir = "NULL"),
          definition = function(ffdir, 
                                ..., 
                                gibbsSamples, 
                                x,
                                spatialData,
                                dB,
                                basisdB,
                                fos) { 

              prep <- .yFinestPrep(spatialData = spatialData,
                                   dB = dB,
                                   fos = fos)

              yFinest <- prep %*% x %*% gibbsSamples$beta +
                         gibbsSamples$xi + 
                         basisdB %*% gibbsSamples$eta

              return( yFinest )
            })

setMethod(f = ".yFinest",
          signature = c(ffdir = "character"),
          definition = function(ffdir, 
                                ..., 
                                gibbsSamples, 
                                x,
                                spatialData,
                                dB,
                                basisdB,
                                fos) { 

              prep <- .yFinestPrep(spatialData = spatialData,
                                   dB = dB,
                                   fos = fos)

              i1 <- NULL
              i2 <- NULL

              yFinest <- ff::ffcolapply(EXPR = prep %*% x %*% gibbsSamples$beta[,i1:i2] +
                                               gibbsSamples$xi[,i1:i2] + 
                                               basisdB %*% gibbsSamples$eta[,i1:i2],
                                        X = gibbsSamples$beta, 
                                        RETURN = TRUE, 
                                        FF_RETURN = TRUE,
                                        RETROW = nrow(x = prep))

              return( yFinest )
})
