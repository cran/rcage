# Generate Plots of Key Results
#
# function is not exported
#
# Generates plots of the average modeled outcome of interest on the finest
#   resolution, the average outcome of interest on the provided clustering,
#   the standard error on the provided clustering, and the cage/DCAGE.
#
#  @param ... Ignored.
#
#  @param dB A Spatial object. The finest resolution. {nB}
#
#  @param yFinest A matrix or ff_matrix object. The modeled outcome 
#    of interest. {nB x nKept}
#
#  @param yOpt A matrix or ff_matrix object. The outcome of interest on
#    the given clustering. {nC x nKept}
#
#  @param cluster A vector object. The cluster ids for each areal unit {nB}
#
#  @param cage A vector object. The CAGE/DCAGE for each areal unit {nB}
#
#  @param dataScale A character object. Function to use if data should be 
#    rescaled before plotting.
#
#  @param palette A character object. The palette of viridis to use.
#
#  @param criterion A character object. CAGE/DCAGE for title.
#
#' @import ggplot2
#' @include summaryInfo.R
#' @importFrom gridExtra  grid.arrange
.mapIt <- function(...,
                   dB,
                   yFinest,
                   yOpt,
                   cluster,
                   cage,
                   dataScale,
                   palette,
                   criterion) {

  # identify unique cluster ids for clustering
  uniqueID <- sort(x = unique(x = cluster))
  nC <- length(x = uniqueID)
  if (nC == 0L) stop("cluster not appropriately defined")

  # number of clusters in optimal clustering
  nClusters <- nrow(x = yOpt)

  # number of Gibbs samples
  nGibbs <- ncol(x = yOpt)

  # number of areal units in finest resolution
  nB <- length(x = cluster)

  # verify that the clustering vector agrees with the finest resolution object
  if (nB != length(x = dB)) {
    stop("cluster is not appropriately dimensioned", call. = FALSE)
  }

  # mean and variance over Gibbs samples for each cluster of optimal
  # list of objects of length {nC}
  sumInfo <- .summaryInfo(y = yOpt, dataScale = dataScale)

  i1 <- NULL
  i2 <- NULL

  if (is(object = yFinest, class2 = "ff_matrix")) {
    if (!is.null(x = dataScale)) {

      yFinest <- ff::ffcolapply(EXPR = do.call(what = dataScale, 
                                               args = list(yFinest[,i1:i2])),
                                X = yFinest, RETURN = TRUE, FF_RETURN = TRUE,
                                CFUN = "ccbind")

      yFinest <- ff::ffrowapply(EXPR = rowMeans(x = yFinest[i1:i2,,drop=FALSE]),
                                X = yFinest, RETURN = TRUE, FF_RETURN = FALSE,
                                CFUN = "c")

    } else {

      yFinest <- ff::ffrowapply(EXPR = rowMeans(x = yFinest[i1:i2,,drop=FALSE]),
                                X = yFinest, RETURN = TRUE, FF_RETURN = FALSE,
                                CFUN = "c")

    }
  } else {

    if (!is.null(x = dataScale)) {

      nyFinest <- ncol(x = yFinest)
      yFinest <- rowMeans(x = matrix(data = do.call(what = dataScale, 
                                                    args = list(yFinest)), 
                                     ncol = nyFinest))

    } else {

      yFinest <- rowMeans(x = yFinest)

    }
  }

  # match cluster values to components that make up the clusters
  orderedMean <- numeric(length = nB)
  orderedSD <- numeric(length = nB)
  ordereDCAGE <- numeric(length = nB)

  for (i in 1L:nC) {
    ind <- cluster == uniqueID[i]
    orderedMean[ind] <- sumInfo$meanY[i]
    orderedSD[ind] <- sumInfo$sd[i]
    ordereDCAGE[ind] <- cage[i]
  }

  # determine if the scale can be simplified
  maxY <- ceiling(x = max(abs(x = orderedMean)))
  idi <- 1L
  while( TRUE ) {
    if (maxY %% 10^idi == abs(x = maxY)) break
    idi <- idi + 1L
  }

  idi <- idi - 1L

  if( idi > 3L ) {
    orderedMean <- orderedMean / {10.0^idi}
    yFinest <- yFinest / {10.0^idi}
    legScale <- paste0("x 10^",idi)
  } else {
    legScale <- ""
  }

  maxY <- ceiling(x = max(abs(x = orderedSD)))

  idi <- 1L
  while( TRUE ) {
    if (maxY %% 10^idi == abs(x = maxY)) break
    idi <- idi + 1L
  }

  idi <- idi - 1L

  if( idi > 2L ) {
    orderedSD <- orderedSD / {10.0^idi}
    legScaleSD <- paste0("x 10^",idi)
  } else {
    legScaleSD <- ""
  }

  # add values to the spatial data data.frame
  df <- data.frame("meanY" = orderedMean,
                   "sqVarY" = orderedSD,
                   "cage" = ordereDCAGE,
                   "yFinest" = yFinest)

  if (is(object = dB, class2 = "SpatialPolygonsDataFrame") ||
      is(object = dB, class2 = "SpatialPointsDataFrame")) {
    dB@data <- df
  } else if (is(object = dB, class2 = "SpatialPolygons")) {
    dB <- sp::SpatialPolygonsDataFrame(Sr = dB, data = df)
  } else if (is(object = dB, class2 = "SpatialPoints")) {
    dB <- sp::SpatialPointsDataFrame(coords = dB, data = df)
  } else {
    stop("inappropriate spatial object")
  }

  # cluster data according to clustering input
  clusteredData <- tryCatch(expr = aggregate(x = dB,
                                             by = list(cluster),
                                             FUN = mean,
                                             dissolve = TRUE,
                                             areaWeighted = FALSE),
                            error = function(e) {
                                      stop("unable to aggregate data to cluster\n", 
                                           e$message, call. = FALSE)
                                    })

  meanY <- NULL
  sqVarY <- NULL

  if (is(object = dB, class2 = "SpatialPolygons")) {

    dB.f <- sf::st_as_sf(x = dB)

    Map1 <- .ggplot2Poly(data = dB.f, 
                         var = yFinest, 
                         legendScale = legScale, 
                         title = "Finest Resolution", 
                         subtitle = "Response", 
                         palette = palette,
                         limits = c(min(dB.f$yFinest), max(dB.f$yFinest)),
                         lwd = 0)

    clusteredData.f <- sf::st_as_sf(x = clusteredData)

    Map2 <- .ggplot2Poly(data = clusteredData.f, 
                         var = meanY, 
                         legendScale = legScale, 
                         title = "Support-Level", 
                         subtitle = "Prediction", 
                         palette = palette,
                         limits = c(min(dB.f$yFinest), max(dB.f$yFinest)),
                         lwd = 0)

    Map3 <- .ggplot2Poly(data = clusteredData.f, 
                         var = sqVarY, 
                         legendScale = legScaleSD, 
                         title = "Support-Level", 
                         subtitle = "Root Prediction Error", 
                         palette = palette,
                         limits = c(min(clusteredData.f$sqVarY), max(clusteredData.f$sqVarY)),
                         lwd = 0)

    Map4 <- .ggplot2Poly(data = clusteredData.f, 
                         var = cage, 
                         legendScale = "", 
                         title = "Support-Level", 
                         subtitle = criterion, 
                         palette = palette,
                         limits = c(min(clusteredData.f$cage), max(clusteredData.f$cage)),
                         lwd = 0)

  } else {

    dB.f <- data.frame(dB)
    dB.f$mmyylloonngg <- sp::coordinates(dB)[,1L]
    dB.f$mmyyllaatt <- sp::coordinates(dB)[,2L]

    Map1 <- .ggplot2Point(data = dB.f, 
                          var = yFinest, 
                          legendScale = legScale, 
                          title = "Finest Resolution", 
                          subtitle = "Response", 
                          palette = palette,
                          limits = c(min(dB.f$yFinest), max(dB.f$yFinest)))

    clusteredData.f <- data.frame(clusteredData)
    clusteredData.f$mmyylloonngg <- sp::coordinates(clusteredData)[,1L]
    clusteredData.f$mmyyllaatt <- sp::coordinates(clusteredData)[,2L]

    Map2 <- .ggplot2Point(data = clusteredData.f, 
                          var = meanY, 
                          legendScale = legScale, 
                          title = "Support-Level", 
                          subtitle = "Prediction", 
                          palette = palette,
                          limits = c(min(dB.f$yFinest), max(dB.f$yFinest)))

    Map3 <- .ggplot2Point(data = clusteredData.f, 
                          var = sqVarY, 
                          legendScale = legScaleSD, 
                          title = "Support-Level", 
                          subtitle = "Root Prediction Error", 
                          palette = palette,
                          limits = c(min(clusteredData.f$sqVarY), max(clusteredData.f$sqVarY)))

    Map4 <- .ggplot2Point(data = clusteredData.f, 
                          var = cage, 
                          legendScale = "", 
                          title = "Support-Level", 
                          subtitle = criterion, 
                          palette = palette,
                          limits = c(min(clusteredData.f$cage), max(clusteredData.f$cage)))
  }

  # Plot results
  plotList <- list(Map1, Map2, Map3, Map4)

  gridExtra::grid.arrange(grobs = plotList, ncol = 2L)
}

#' @importFrom rlang enquo
#' @import sf
#
.ggplot2Poly <- function(data, 
                         var, 
                         legendScale, 
                         title, 
                         subtitle, 
                         palette, 
                         limits, 
                         lwd) {

  long <- NULL
  lat <- NULL
  group <- NULL
  quo_var <- rlang::enquo(var)

  Map <- ggplot2::ggplot(data, ggplot2::aes(x = long, 
                                            y = lat,  
                                            group = group,  
                                            fill = !! quo_var))

  Map <- Map + ggplot2::geom_sf(lwd = lwd)

  Map <- Map + ggplot2::coord_sf()

  Map <- Map + ggplot2::labs(fill = legendScale)

  Map <- Map + 
         ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                        axis.title.y = ggplot2::element_blank(),
                        axis.text.x = ggplot2::element_blank(),
                        axis.text.y = ggplot2::element_blank(),
                        axis.ticks.x = ggplot2::element_blank(),
                        axis.ticks.y = ggplot2::element_blank(),
                        legend.position = "bottom",
                        legend.title = ggplot2::element_text(),
                        plot.title = ggplot2::element_text(size = 10))

  Map <- Map + ggplot2::ggtitle(label = title, subtitle = subtitle)

  Map <- Map + ggplot2::scale_fill_viridis_c(limits = limits, option = palette)

  return( Map )
}

#' @importFrom rlang enquo

.ggplot2Point <- function(data, 
                          var, 
                          legendScale, 
                          title, 
                          subtitle, 
                          palette, 
                          limits) {

  mmyylloonngg <- NULL
  mmyyllaatt <- NULL
  quo_var <- rlang::enquo(var)

  Map <- ggplot2::ggplot(data, ggplot2::aes(x = mmyylloonngg, 
                                            y = mmyyllaatt,  
                                            color = !! quo_var)) + 
         ggplot2::geom_point() +  
         ggplot2::coord_equal() +  
         ggplot2::labs(fill = legendScale) +  
         ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                        axis.title.y = ggplot2::element_blank(),
                        axis.text.x = ggplot2::element_blank(),
                        axis.text.y = ggplot2::element_blank(),
                        axis.ticks.x = ggplot2::element_blank(),
                        axis.ticks.y = ggplot2::element_blank(),
                        legend.position = "bottom",
                        legend.title = ggplot2::element_text(),
                        plot.title = ggplot2::element_text(size = 10)) + 
         ggplot2::ggtitle(label = title, subtitle = subtitle)

  Map <- Map + ggplot2::scale_color_viridis_c(limits = limits, option = palette)

  return( Map )
}
