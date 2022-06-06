#' Generate Plots Summarizing Main Results
#'
#' Function generates plots showing the finest resolution, the support-level
#'  predictions, the support-level root prediction error, and the support-level
#'  CAGE/DCAGE.
#'
#' @param x An rcage object. The value returned by optRegion() or region().
#'
#' @param ... ignored.
#'
#' @param dB A SpatialPolygons object or SpatialPoints object. The areal
#'    units of the finest resolution.
#'
#' @param dataScale A function or NULL. If not NULL, the function to be 
#'   applied to the outcome interest to adjust the scaling in plots.
#'
#' @param palette A character. The palette preference for plotting. The 
#'   palette is assumed to be viridis palette with possible values
#'   'viridis', 'magma', 'plasma', 'inferno', 'cividis'
#'
#' @method plot rcage
#' @export
#'
#' @returns No value object returned; called to generate plots.
#'
#' @examples
#' 
#' # create 5x5 square 
#' 
#' poly <- raster::rasterToPolygons(raster::raster(nrows = 5, ncols = 5,
#'                                                 xmn = -1.25, xmx = 1.25,
#'                                                 ymn = -1.25, ymx = 1.25,
#'                                                 res = 0.5,
#'                                                 crs = "+proj=longlat +datum=WGS84"))
#' 
#' df <- data.frame("x" = stats::rnorm(n = 25))
#' 
#' dt <- sp::SpatialPolygonsDataFrame(poly, df)
#' 
#' knots <- cbind(c(-0.75, 0.0, 0.75, -0.75, 0.0, 0.75, -0.75, 0.0, 0.75),
#'                c(-0.75,-0.75, -0.75, 0.0, 0.0, 0.0, 0.75, 0.75, 0.75))
#'
#' res <- optRegion(spatialData = dt,
#'                  response = "x",
#'                  sigmavar = rep(1, 25),
#'                  gL = 5,
#'                  gU = 7,
#'                  nGibbs = 50L,
#'                  nBurn = 10L,
#'                  nThin = 1L,
#'                  nw = 2000L,
#'                  knots = knots)
#'
#'  plot(x = res, dB = dt)

plot.rcage <- function(x, ..., 
                       dB,
                       dataScale = NULL,
                       palette = "plasma") {

  if (is.character(x = x$yFinest)) {
    yFinest <- NULL
    ffload(file = x$yFinest, overwrite = TRUE)
    x$yFinest <- yFinest
  }

  if (is.character(x = x$yOpt)) {
    yOpt <- NULL
    ffload(file = x$yOpt, overwrite = TRUE)
    x$yOpt <- yOpt
  }

  if (length(x = dB) != nrow(x = x$yFinest)) stop("dim of dB incorrect")

  .mapIt(dB = dB,
         yFinest = x$yFinest,
         yOpt = x$yOpt,
         cluster = x$cluster$cluster,
         cage = x$CAGETrack,
         dataScale = dataScale,
         palette = palette,
         criterion = x$criterion)

}
