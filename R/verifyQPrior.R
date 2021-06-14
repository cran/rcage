# @param ... ignored.
#
# @param qPrior A character object or NULL. The name of the prior to be
#   used. Must be one of 'mi' or 'wish' indicating the Moran I prior or the
#   Inverse Wishart. If NULL, the MI prior is assumed.
#
# @return A default QPriorObj_MI or QPriorObj_Wishart object.
#
#' @include QPriorObj.R QPriorObj_MI.R QPriorObj_Wishart.R
.verifyQPrior <- function(..., qPrior) { 

  if (is.null(x = qPrior) || grepl(pattern = "mi", x = tolower(x = qPrior))) {

    return( new(Class = "MI") ) 

  } else if (grepl(pattern = "wish", x = tolower(x = qPrior))) {

    return( new(Class = "Wishart") )
    
  } else {

    stop("unrecognized Qprior", call. = FALSE)

  }
}
