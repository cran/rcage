#  Calculate P_W Lambda_W ^{-1/2} - Cholesky Decomposition of W^{-1}
#
#  function is not exported
#
#  @param ... Ignored. Included only to require names inputs.
#
#  @param spatialData A SpatialPoints object, SpatialPolygons object, or 
#    list of said objects. The source support.
#
#  @param dB A SpatialPolygons object, integer, or NULL. The finest resolution.
#
#  @param gbfObj A GBFObj object. Information about basis function.
#
#  @param nw A numeric object or NULL. MC replicates.
#
#  @param verify A logical object. If TRUE, verify non-singular matrix.
#
#  @return a list containing
#     matrix : Pw Lambdaw^{-1/2}  { r x r }, 
#     w : final bandwidth
#
#' @include wMatrix.R
#
.cholesky_W <- function(...,
                        spatialData, 
                        dB,
                        gbfObj,  
                        nw,
                        verify) {

  message("\tcalculating Cholesky square root of W^{-1}")

  success <- FALSE
  cnt <- 0L
  rr <- gbfObj@args$w
  reset <- FALSE

  while( cnt < 100L ) {

    # calculates W Matrix {r x r}
    # only print information to the screen on the first iteration
    wMatrix <- .wMatrix(spatialData = spatialData, 
                        nw = nw,  
                        dB = dB,
                        gbfObj = gbfObj,
                        verbose = {cnt == 0L})

    # W^{-1} {r x r}
    # attempt to invert W
    wInverse <- tryCatch(expr = solve(a = wMatrix),
                         error = function(e) { return(e$message) })

    # if inverse was successful on first try, use provided w
    if (is.matrix(x = wInverse) && cnt == 0L) break

    # if inverse was twice successful use current value of w
    if (is.matrix(x = wInverse) && success) break

    # if inverse was finally successful set success and increase w
    # goal is to have w slightly higher than "minimum"
    if (is.matrix(x = wInverse)) {
      success <- TRUE
    } else {
      success <- FALSE
    }

    # if not modifying/verifying w, stop with error msg given by solve()
    if (!verify) stop(wInverse, call. = FALSE)

    # cap the number of attempts at 100
    if (cnt < 99L) {
      if (any(colSums(x = abs(x = wMatrix)) <= 1e-8)) {
        # if any columns are all 0, increase basis function bandwidth and
        # try again
        gbfObj@args$w <- 1.01*gbfObj@args$w
        rr <- gbfObj@args$w
        reset <- TRUE
      } else if (success) {
        # if first successful inverse, increase basis function bandwidth only
        # slightly and try again
        gbfObj@args$w <- 1.001*gbfObj@args$w
        rr <- gbfObj@args$w
      } else {
        # if not successful on second try, revert to successful value
        # and exit
        gbfObj@args$w <- gbfObj@args$w / 1.001
        rr <- gbfObj@args$w
        break
      }
      cnt <- cnt + 1L
    } else {
      stop("unable to identify an appropriate bandwidth\n", 
           wInverse, call. = FALSE)
    }
  }

  if (reset) {
    message("\tunable to invert W at original bandwidth; increased to ", rr)
  }

  # Cholesky decomposition of W^{-1} { r x r }
  cholesky_W <- tryCatch(expr = chol(x = wInverse),
                         error = function(e) {
                                   stop("cholesky decomposition of W^{-1} failed\n",
                                        e$message, call. = FALSE)
                                 })

  return( list("matrix" = t(x = cholesky_W), "w" = gbfObj@args$w) )

}
