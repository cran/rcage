# Class to define QPrior functions
#
# Class and methods are not exported
#
setClass(Class = "QPriorObj", 
         slots = c("Q" = "matrix", 
                   "deno" = "matrix", 
                   "mu" = "ANY",
                   "lambda" = "numeric"))

setMethod(f = "initialize",
          signature = "QPriorObj",
          definition = function(.Object, ...){
              .Object@Q <- matrix(data = 0.0, nrow = 0L, ncol = 0L)
              .Object@deno <- matrix(data = 0.0, nrow = 0L, ncol = 0L)
              .Object@mu <- matrix(data = 0.0, nrow = 0L, ncol = 0L)
              .Object@lambda <- numeric(length = 0L)
               
              return( .Object )
            })

setGeneric(name = ".qInv",
           def = function(qObj, ...) { standardGeneric(".qInv") })

setMethod(f = ".qInv",
          signature = c(qObj = "ANY"),
          definition = function(qObj, ...) { stop("not allowed") })


setGeneric(name = ".gibbsQ",
           def = function(qObj, ...) { standardGeneric(".gibbsQ") })

setMethod(f = ".gibbsQ",
          signature = c(qObj = "ANY"),
          definition = function(qObj, ...) { stop("not allowed") })

setGeneric(name = ".metroQ",
           def = function(qObj, ...) { standardGeneric(".metroQ") })

setMethod(f = ".metroQ",
          signature = c(qObj = "ANY"),
          definition = function(qObj, ...) { stop("not allowed") })

.gibbsStep <- function(..., qObj, denoEta, SpinvV) {

  msg <- "error/warning in Cholesky decomposition for eta distribution\n"

  # Cholesky decomposition of {Qinv / lambda + SpinvV*S}^{-1}
  qObj@deno <- tryCatch(expr = chol(x = denoEta), 
                        condition = function(e){
                                      stop(msg, e$message, call. = FALSE)
                                    })

  qObj@mu <- denoEta %*% SpinvV

  return( qObj )

}
