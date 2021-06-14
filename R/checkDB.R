# S4 methods to simply the verification of input dB
#
# Allowed combinations and their returned values
#
# spatialData	dB	finest			dB	cage
# input		input   is taken as input       out	out
# ------------------------------------------------------------------
# pts		NULL	spatialData		NULL	TRUE
# pts		poly	dB			dB	TRUE/FALSE
# poly		NULL	spatialData		NULL	FALSE
# poly		poly	dB			dB	FALSE
# list		NULL	error
# list		int (i)	m poly n pts
#			n>1 error
#			m=0 n=1 error
#			m=1 n=0 error
#			all others list[[i]]	i	TRUE/FALSE
# list		poly	m poly n pts
#			n>1 error
#			m=0 n=1 error
#			m=1 n=0 error
#			all other		dB	TRUE/FALSE


#' @import sp
#' @import methods
setGeneric(name = ".checkDB",
           def = function(spatialData, dB, ...) { standardGeneric(".checkDB") })

# anything not explicitly allowed is forbidden
setMethod(f = ".checkDB",
          signature = c(spatialData = "ANY",
                        dB = "ANY"),
          definition = function(spatialData, dB, ...) { stop("not allowed") })

# SpatialPoints and NULL indicates Ds is the finest resolution.
# CAGE is the only option
setMethod(f = ".checkDB",
          signature = c(spatialData = "SpatialPoints",
                        dB = "NULL"),
          definition = function(spatialData, dB, ...) {

              message("source data are SpatialPoints (D_S)\n",
                      "D_B = D_S\n", "CAGE")

              return( list("dB" = NULL, "cage" = TRUE) )
            })

# SpatialPoints w/ a SpatialPolygons indicates that dB is the finest resolution.
# Can be either CAGE or DCAGE depending on user input.
#' @importFrom rgeos gContains
setMethod(f = ".checkDB",
          signature = c(spatialData = "SpatialPoints",
                        dB = "SpatialPolygons"),
          definition = function(spatialData, dB, ..., cage) {

              # if provided, the finest resolution must contain all source
              # data
              if (!rgeos::gContains(spgeom1 = dB, spgeom2 = spatialData)) {
                stop("dB does not fully contain spatialData", call. = FALSE)
              }

              message("source data are SpatialPoints\n",
                      "D_B is SpatialPolygons and provided as input\n",
                      ifelse(test = cage, yes = "CAGE", no = "DCAGE"))

              # limit dB to only those components of dB for which point
              # data are provided
              ovr <- rgeos::gContains(spgeom1 = dB, 
                                      spgeom2 = spatialData,  
                                      byid = TRUE)

              if (any(colSums(x = ovr) == 0L)) {
                dB <- dB[colSums(x = ovr) > 0L,,drop = FALSE]
                message("removed ", sum(colSums(x = ovr) == 0L), 
                        " elements from D_B as they do not overlap the data")
              }

              return( list("dB" = dB, "cage" = cage) )
            })

# SpatialPolygons w/ NULL indicates that D_A is the finest resolution. 
# DCAGE is the only option.
setMethod(f = ".checkDB",
          signature = c(spatialData = "SpatialPolygons",
                        dB = "NULL"),
          definition = function(spatialData, dB, ...) {

              message("source data are SpatialPolygons (D_A)\n",
                      "D_B = D_A\n", "DCAGE")

              return( list("dB" = NULL, "cage" = FALSE) )
            })

# SpatialPolygons w/ SpatialPolygons indicates that dB is the finest resolution.
# DCAGE is the only option.
setMethod(f = ".checkDB",
          signature = c(spatialData = "SpatialPolygons",
                        dB = "SpatialPolygons"),
          definition = function(spatialData, dB, ...) {

              # if provided, the finest resolution must contain all source
              # data
              if (!rgeos::gContains(spgeom1 = dB, spgeom2 = spatialData)) {
                stop("dB does not fully contain spatialData", call. = FALSE)
              }

              message("source data are SpatialPolygons (D_A)\n",
                      "D_B is SpatialPolygons and provided as input\n", 
                      "DCAGE")

              # limit dB to only those components of dB for which poly
              # data are provided
              ovr <- rgeos::gContains(spgeom1 = dB, 
                                      spgeom2 = spatialData,  
                                      byid = TRUE)
              if (any(colSums(x = ovr) == 0L)) {
                dB <- dB[colSums(x = ovr) > 0L,drop = FALSE]
                message("removed ", sum(colSums(x = ovr) == 0L), 
                        " elements from D_B as they do not overlap the data")
              }

              return( list("dB" = dB, "cage" = FALSE) )
            })

# list and NULL is not allowed.  
setMethod(f = ".checkDB",
          signature = c(spatialData = "list",
                        dB = "NULL"),
          definition = function(spatialData, dB, ..., cage) {
              stop("if spatialData is a list, dB must be either\n",
                   "\tan integer specifying the element to be used as the ",
                     "finest resolution, or\n",
                   "\ta SpatialPolygons object defining the finest resolution", 
                   call. = FALSE)
            })

# list w/ SpatialPolygons indicates that dB is the finest resolution
setMethod(f = ".checkDB",
          signature = c(spatialData = "list",
                        dB = "SpatialPolygons"),
          definition = function(spatialData, dB, ..., cage) {

              if (length(x = spatialData) == 1L) {
                stop("spatialData cannot be a list of length 1", call. = FALSE)
              }

              # count number of SpatialPoints and SpatialPolygons objects
              npt <- 0L
              nply <- 0L

              for (i in 1L:length(x = spatialData)) {

                # ensure the dB contains all source data
                if (!rgeos::gContains(spgeom1 = dB, 
                                      spgeom2 = spatialData[[ i ]])) {
                  stop("dB does not fully contain spatialData", 
                       call. = FALSE)
                }

                if (is(object = spatialData[[ i ]], class2 = "SpatialPoints")) {
                  npt <- npt + 1L
                } else if (is(object = spatialData[[ i ]], 
                              class2 = "SpatialPolygons")) {
                  nply <- nply + 1L
                } else {
                  stop("unrecognized object provided in spatialData", 
                       call. = FALSE)
                }
              }

              # if no SpatialPoints data provided, must be DCAGE
              if (npt == 0L) cage <- FALSE

              if (npt > 1L) {
                stop("multiple SpatialPoints objects provided as support ",
                     "combine before calling this method", call. = FALSE)
              }

              if (npt == 1L) {

                message("source data are 1 SpatialPoints (D_S) and ", nply,
                        " SpatialPolygons (D_A)\n",
                        "D_B is SpatialPolygons provided as input\n", 
                        ifelse(test = cage, yes = "CAGE", no = "DCAGE"))

              } else {

                message("source data are ", nply, " SpatialPolygons (D_A)\n",
                        "D_B is SpatialPolygons provided as input\n", 
                        "DCAGE")

              }

              # limit dB to only those components of dB for which 
              # data are provided
              keep <- numeric(length = length(x = dB))
              for (i in 1L:length(x = spatialData)) {
                ovr <- rgeos::gContains(spgeom1 = dB, 
                                        spgeom2 = spatialData[[ i ]],  
                                        byid = TRUE)
                keep <- keep + colSums(x = ovr)
              }
              if (any(keep == 0L)) {
                dB <- dB[keep > 0,,drop = FALSE]
                message("removed ", sum(keep == 0L), 
                        " elements from D_B as they do not overlap the data")
              }

              return( list("dB" = dB, "cage" = cage) )

            })

# list w/ integer indicates that element dB is the finest resolution
setMethod(f = ".checkDB",
          signature = c(spatialData = "list",
                        dB = "numeric"),
          definition = function(spatialData, dB, ..., cage) {

              if (length(x = spatialData) == 1L) {
                stop("spatialData cannot be a list of length 1", call. = FALSE)
              }

              if (dB <= 0L || dB > length(x = spatialData)) {
                stop("inappropriate value provided for dB", call. = FALSE)
              }

              if (is(object = spatialData[[ dB ]], class2 = "SpatialPoints")) {
                stop("dB cannot be SpatialPoints", call. = FALSE)
              }

              # count number of SpatialPoints and SpatialPolygons objects
              npt <- 0L
              nply <- 0L

              for (i in 1L:length(x = spatialData)) {

                if (i != dB) {
                  # ensure the dB contains all source data
                  if (!rgeos::gContains(spgeom1 = spatialData[[ dB ]], 
                                        spgeom2 = spatialData[[ i ]])) {
                    stop("dB does not fully contain spatialData", 
                         call. = FALSE)
                  }
                }

                if (is(object = spatialData[[ i ]], class2 = "SpatialPoints")) {
                  npt <- npt + 1L
                } else if (is(object = spatialData[[ i ]], 
                              class2 = "SpatialPolygons")) {
                  nply <- nply + 1L
                }
              }

              # if no SpatialPoints data provided, must be DCAGE
              if (npt == 0L) cage <- FALSE

              if (npt > 1L) {
                stop("multiple SpatialPoints objects provided as support ",
                     "combine before calling this method", call. = FALSE)
              }

              if (npt == 1L) {

                message("source data are 1 SpatialPoints (D_S) and ", nply,
                        " SpatialPolygons (D_A)\n",
                        "D_B is element ", dB, " of spatialData\n", 
                        ifelse(test = cage, yes = "CAGE", no = "DCAGE"))

              } else {

                message("source data are ", nply, " SpatialPolygons (D_A)\n",
                        "D_B is element ", dB, " of spatialData\n", 
                        "DCAGE")

              }

              return( list("dB" = dB, "cage" = cage) )

            })
