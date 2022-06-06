#  Mean and Variance Over Gibbs Samples
#
#  method is not exported
#
#  Calculates mean an standard error for plot method
#
#  @param y A matrix or ff_matrix object. The modeled outcome of interest.
#
#  @param dataScale A character or function object or NULL. The 
#    function to use to scale data prior to plotting.
#
#  @param ... ignored.
#
#  @return A list containing
#   \item{meanY} vector of mean of y {nC} 
#   \item{sd} vector of standard deviation of y {nC} 
setGeneric(name = ".summaryInfo", 
           def = function(y, dataScale, ...) { 
                     standardGeneric(".summaryInfo") 
                   })

# default method results in error
setMethod(f = ".summaryInfo",
          signature = c(y = "ANY",
                        dataScale = "ANY"),
          definition = function(y, dataScale, ...) { stop("not allowed") })

# call scaling function on data then call method with no scaling
setMethod(f = ".summaryInfo",
          signature = c(y = "matrix",
                        dataScale = "ANY"),
          definition = function(y, dataScale, ...) { 
              y <- do.call(what = dataScale, args = list(y))
              return( .summaryInfo(y = y, dataScale = NULL) ) 
            })

# summary statistics
#' @importFrom stats var
setMethod(f = ".summaryInfo",
          signature = c(y = "matrix",
                        dataScale = "NULL"),
          definition = function(y, dataScale, ...) { 
              meanY <- rowMeans(x = y)
              varY  <- apply(X = y, MARGIN = 1L, FUN = stats::var)
              sqVarY <- sqrt(x = varY)
              return( list("meanY" = meanY, "sd" = sqVarY) )
            })

# call scaling function on data then call method with no scaling
#' @import ff
setMethod(f = ".summaryInfo",
          signature = c(y = "ff_matrix",
                        dataScale = "ANY"),
          definition = function(y, dataScale, ...) { 
              # only to appease cran testing
              i1 <- NULL
              i2 <- NULL

              yS <- ff::ffrowapply(do.call(what = dataScale, 
                                           args = list(y[i1:i2,])), 
                                   X = y,
                                   CFUN = 'crbind',
                                   RETURN = TRUE,
                                   FF_RETURN = TRUE)

              return( .summaryInfo(y = yS, dataScale = NULL) ) 

            })

# summary statistics
#' @import ff
#' @importFrom stats var
setMethod(f = ".summaryInfo",
          signature = c(y = "ff_matrix",
                        dataScale = "NULL"),
          definition = function(y, dataScale, ...) { 

              # only to appease cran testing
              i1 <- NULL
              i2 <- NULL

              meanY <- ff::ffrowapply(rowMeans(x = y[i1:i2,]), 
                                      X = y,
                                      CFUN = 'c',
                                      RETURN = TRUE,
                                      FROM = "i1",
                                      TO = "i2",
                                      FF_RETURN = FALSE)

              varY <- ff::ffrowapply(apply(X = y[i1:i2,], 
                                           MARGIN = 1L, 
                                           FUN = stats::var), 
                                     X = y,
                                     CFUN = 'c',
                                     RETURN = TRUE,
                                     FROM = "i1",
                                     TO = "i2",
                                     FF_RETURN = FALSE)

              sdY <- sqrt(x = drop(x = varY))

              return( list("meanY" = drop(x = meanY), "sd" = sdY) )

            })
