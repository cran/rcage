#' @include QPriorObj.R

#  Information required for Inverse Wishart Prior
#
#  @slot scale A "matrix" object. The scale matrix for the distribution.
#
setClass(Class = "Wishart", 
         slots = c("scale" = "matrix"),
         contains = c("QPriorObj"))

.validityWishart <- function(object) {

  if (length(x = object@lambda) != ncol(x = object@Q)) {
    return("dim error lambda")
  }

  if (nrow(x = object@scale) != ncol(x = object@Q)) {
    return("dim error scale")
  }

  if (ncol(x = object@scale) != ncol(x = object@Q)) {
    return("dim error scale")
  }

  if (any(object@scale < 0.0)) return("scale cannot be negative")

  return( TRUE )

}

setValidity(Class = "Wishart", method = .validityWishart)

setMethod(f = "initialize",
          signature = "Wishart",
          definition = function(.Object, ...){

              .Object <- callNextMethod()
              .Object@scale <- matrix(data = 0.0, nrow = 0L, ncol = 0L)

              validObject(.Object)

              return( .Object )
            })

# @param basisdB is a matrix of the phiOC {nB x r}
# @param lambda is a numeric vector of the eigenvalues
# @param scale a numeric vector
setMethod(f = ".qInv",
          signature = c(qObj = "Wishart"),
          definition = function(qObj, ..., basisdB, lambda, scale) { 

              message("initializing Wishart prior")

              if (!is.numeric(x = scale)) stop("wishartScale must be numeric", 
                                               call. = FALSE)

              r <- ncol(x = basisdB)

              qObj@Q <- diag(x = 1.0, nrow = r, ncol = r)

              if (length(x = scale) != 1L && length(x = scale) != r) {
                stop("wishartScale must be of length 1 or r", call. = FALSE)
              }
              qObj@scale <- diag(x = scale, nrow = r, ncol = r)


              if (length(x = lambda) != 1L && length(x = lambda) != r) {
                stop("lambda must be of length 1 or r", call. = FALSE)
              }

              if (length(x = lambda) == 1L) lambda <- rep(x = lambda, times = r)

              qObj@lambda <- lambda

              return( qObj )

            })

#' @importFrom pracma pinv
setMethod(f = ".gibbsQ",
          signature = c(qObj = "Wishart"),
          definition = function(qObj, ..., SpinvVS, SpinvV) {

              r <- ncol(x = qObj@Q)

              denoEta <- as.matrix(x = SpinvVS) + 
                         qObj@Q %*% {{1.0 / qObj@lambda} * t(x = qObj@Q)}

              msg <- "unable to invert matrix for eta distribution\n"

              denoEta <- tryCatch(expr = pracma::pinv(A = as.matrix(x = denoEta)), 
                                  error = function(e){
                                            stop(msg, e$message, call. = FALSE)
                                          })

              qObj <- .gibbsStep(qObj = qObj, 
                                 denoEta = denoEta, 
                                 SpinvV = SpinvV)

              return( qObj )
            })


#' @importFrom LaplacesDemon rinvwishart
setMethod(f = ".metroQ",
          signature = c(qObj = "Wishart"),
          definition = function(qObj, ..., eta) {

              r <- length(x = eta)

              scale <- qObj@scale + tcrossprod(x = as.matrix(x = eta))

              matsim <- LaplacesDemon::rinvwishart(nu = {r + 1L}, S =  scale)

              if (any(is.na(x = matsim)) || 
                  any(is.nan(x = matsim)) || 
                  any(is.infinite(x = matsim))) {
                stop('Wishart contains invalid values (Inf, NA, or NaN)')
              }

              temp <- tryCatch(expr = eigen(x = matsim),
                               error = function(e){
                                         stop("eigen decomposition of prior failed", 
                                              e$message, call. = FALSE)
                                       })
              qObj@Q <- temp$vectors

              return( qObj )
 })
