#' @include QPriorObj.R
#
#  Methods are not exported
#
#  Information required to calculate MI prior
#
# @slots set A logical object. TRUE indicates that matrix has been set.
#   and does not need to be recalculated
#
setClass(Class = "MI", 
         slots = c("set" = "logical"), 
         contains = c("QPriorObj"))

setMethod(f = "initialize",
          signature = "MI",
          definition = function(.Object, ...){

              .Object <- callNextMethod()
              .Object@set <- FALSE

              return( .Object )
            })


#' @include MIpriorInv.R
# @param basisdB is a matrix of the phiOC {nB x r}
# @param lambda is a numeric vector of the eigenvalues
# @param spatialData the source support
# @param dB is either NULL, integer, or SpatialPolygons
# @param dMax is NULL or a numeric indicating the maximum distance to be
#   considered a neighbor (used only when dB is SpatialPoints)
setMethod(f = ".qInv",
          signature = c(qObj = "MI"),
          definition = function(qObj, 
                                ...,
                                basisdB,
                                lambda,
                                spatialData,
                                dB, 
                                dMax) {

              message("initializing MI prior")

              # this is actually the inverse of Q

              if (is.null(x = dB)) {
                # spatialData is the finest resolution
                qObj@Q <- .miPrior(psi = basisdB, dB = spatialData, dMax = dMax)
              } else if (is.numeric(x = dB)) {
                # spatialData is a list and element dB is the finest resolution
                qObj@Q <- .miPrior(psi = basisdB,  
                                   dB = spatialData[[ dB ]], 
                                   dMax = dMax)
              } else {
                # finest resolution provided as a SpatialPolygons object
                qObj@Q <- .miPrior(psi = basisdB, dB = dB, dMax = dMax)
              }

              qObj@lambda <- lambda

              return( qObj )

            })

#' @importFrom pracma pinv
setMethod(f = ".gibbsQ",
          signature = c("qObj" = "MI"),
          definition = function(qObj, ..., SpinvVS, SpinvV) {

              # calculation has already been completed, return
              if (qObj@set) return( qObj )

              msg <- "unable to invert matrix for eta distribution\n"

              # as.matrix is necessary because SpinvVS might be Matrix object
              # {Qinv / lambda + SpinvV*S}^{-1}
              denoEta <- tryCatch(expr = pracma::pinv(A = as.matrix(x = SpinvVS) + 
                                                          qObj@Q / qObj@lambda), 
                                  condition = function(e){
                                                stop(msg, e$message, 
                                                     call. = FALSE)
                                               })

              qObj <- .gibbsStep(qObj = qObj, 
                                 denoEta = denoEta, 
                                 SpinvV = SpinvV)

              qObj@set <- TRUE

              return( qObj )
            })

# metroq step is not performed for MI prior
setMethod(f = ".metroQ",
          signature = c("qObj" = "MI"),
          definition = function(qObj, ...) { return( qObj ) })
