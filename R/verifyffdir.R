# Verify that ff Directory is Appropriately Specified
#
# method not exported
#
# Checks validity of path of ffdir provided. Ensure sufficient memory otherwise.
#
# @param ffdir A character object or NULL. The path to store data.
#
# @param ... ignored.
#
# @param nR An integer object. The maximum number of rows anticipated.
#
# @param nKept An integer object. The maximum number of columns anticipated.
#
# @return NULL
#
# @name verifyffdir
#
#' @import ff
setGeneric(name = ".verifyffdir",
           def = function(ffdir, ...) { standardGeneric(".verifyffdir") })

# default method results in error
setMethod(f = ".verifyffdir",
          signature = c(ffdir = "ANY"),
          definition = function(ffdir, ...) { stop("not allowed") })

setMethod(f = ".verifyffdir",
          signature = c(ffdir = "NULL"),
          definition = function(ffdir, ..., nKept, nR) {

              msg <- "insufficient memory available - define ffdir\n"

              testnThin <- tryCatch(expr = matrix(data = 1.0, 
                                                  nrow = nR,  
                                                  ncol = nKept),
                                    error = function(e) {
                                              stop(msg, e$message, call. = FALSE)
                                            })
              rm(testnThin)
              gc()
           
              return()

            })

setMethod(f = ".verifyffdir",
          signature = c(ffdir = "character"),
          definition = function(ffdir, ...) {
              if (!dir.exists(paths = ffdir)) {
                stop("ffdir not in path", call. = FALSE)
              }
              return()
            })
