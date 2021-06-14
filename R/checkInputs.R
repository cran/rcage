# input cMethod is a character indicating if the hierarchical cluster method
# or the k-means cluster method is requested by the user
.check_cMethod <- function(..., cMethod) {

  # cMethod must be a character
  if (!is.character(x = cMethod)) {
    stop("cMethod must be a character", call. = FALSE)
  }

  # currently limited to 'kmeans' or 'hier'
  cMethod <- tolower(x = cMethod)
  if (!{cMethod %in% c('kmeans', 'hier')}) {
    stop("cMethod must be one of 'kmeans' or 'hier'", call. = FALSE)
  }

  return( cMethod )

}


# ensures that the number of cores is appropriately specified
# and if parallel methods requested, initiates the cluster
# returns NULL of parallel methods are not requested, returns
# the cluster if they are requested.
#' @importFrom parallel makeCluster
.check_nCore <- function(..., nCore, parallelLog) {

  if (is.null(x = nCore)) return( NULL )

  if (!is.numeric(x = nCore)) stop("nCore must be numeric", call. = FALSE)

  nCore <- as.integer(x = nCore)

  if (nCore <= 1L) return( NULL )

  # create cluster include parallelLog specification

  args <- list()
  args[[ "spec" ]] <- nCore
  args[[ "outfile" ]] <- parallelLog

  localCluster <- tryCatch(expr = do.call(what = parallel::makeCluster,
                                          args = args),
                           error = function(e){
                                     stop("unable to create cluster\n",
                                          e$message, call. = FALSE)
                                   })

  # load the ff package across the cluster
  parallel::clusterEvalQ(cl = localCluster, 
                         suppressPackageStartupMessages(expr = library(ff)))

  return( localCluster )
}
