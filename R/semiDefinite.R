# Approximate the best Hermitian positive semi-definite to A in the 2-norm
#
# function is not exported
#
# @param A A matrix object. The matrix for which approximation is sought
#
# @param ... ignored.
#
# @return A matrix object. The best Hermitian PSD approx to A
#
#' @importFrom pracma eps
.semiDefinite <- function(A, ..., verbose = TRUE){

  if (verbose) message("approximating positive semi-definite")

  # symmetrize A
  A <- {A + t(x = A)} / 2.0

  # compute the symmetric polar factor of A.
  svdResult <- tryCatch(expr = svd(x = A),
                        error = function(e) {
                                  stop("unable to obtain PSD approximation\n", 
                                       e$message, call. = FALSE)
                                })

  H <- svdResult$v %*% diag(x = svdResult$d) %*% t(x = svdResult$v)

  # the best Hermitian positive semi-definite approximation
  # to A in the 2-norm
  # http://eprints.ma.man.ac.uk/694/01/covered/MIMS_ep2007_9.pdf
  A <- {A + H} / 2.0

  # ensure symmetry
  A <- {A + t(x = A)} / 2.0

  # test that A is in fact PD. if it is not so, likely due to
  # floating point garbage. Remove small values from diag until a
  # PD is obtained. Stop if not successful after 10 attempts
  k <- 0L
  fixed <- FALSE
  added <- 0.0

  for (p in 1L:10L) {

    res <- tryCatch(expr = chol(x = A),
                    condition = function(e){ FALSE })

    if (is.logical(x = res)) {
      vals <- eigen(x = A, only.values = TRUE)$values

      minEig <- min(Re(z = vals))

      rmItm <- {-minEig * k * k + 
                pracma::eps(x = minEig)} * diag(nrow = nrow(x = A))
      added <- added + rmItm[1L,1L]

      if (added > 1e-4) {
        stop("modification required to make matrix positive definite is too large",
             call. = FALSE)
      }

      A <- A + rmItm

      k <- k + 1L

    } else {

      fixed <- TRUE

      break

    }

  }

  if (!fixed) stop("unable to ensure W is positive definite", call. = FALSE)

  if (k > 0L) {
    message("\t", added, 
            " added to diagonal to ensure positive definiteness")
  }

  return( A )
}
