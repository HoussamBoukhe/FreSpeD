#' @name col_scale
#' @title create a color scale legend for a plot
#' @description  create a color scale legend for a plot 
#' @param z numeric vector, the input vector
#' 
#' @return a color scale legend
col_scale<-function (z, zlim, col = heat.colors(12), breaks, horiz = TRUE, 
                     ylim = NULL, xlim = NULL, ...) 
{
  if (!missing(breaks)) {
    if (length(breaks) != (length(col) + 1)) {
      stop("must have one more break than colors")
    }
  }
  if (missing(breaks) & !missing(zlim)) {
    breaks <- seq(zlim[1], zlim[2], length.out = (length(col) + 
                                                    1))
  }
  if (missing(breaks) & missing(zlim)) {
    zlim <- range(z, na.rm = TRUE)
    zlim[2] <- zlim[2] + c(zlim[2] - zlim[1]) * (0.001)
    zlim[1] <- zlim[1] - c(zlim[2] - zlim[1]) * (0.001)
    breaks <- seq(zlim[1], zlim[2], length.out = (length(col) + 
                                                    1))
  }
  poly <- vector(mode = "list", length(col))
  for (i in seq(poly)) {
    poly[[i]] <- c(breaks[i], breaks[i + 1], breaks[i + 1], 
                   breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if (horiz) {
    YLIM <- c(0, 1)
    XLIM <- range(breaks)
  }
  if (!horiz) {
    YLIM <- range(breaks)
    XLIM <- c(0, 1)
  }
  if (missing(xlim)) 
    xlim = XLIM
  if (missing(ylim)) 
    ylim = YLIM
  plot(1, 1, t = "n", ylim = ylim, xlim = xlim, xaxt = xaxt, 
       yaxt = yaxt, xaxs = "i", yaxs = "i", ...)
  for (i in seq(poly)) {
    if (horiz) {
      polygon(poly[[i]], c(0, 0, 1, 1), col = col[i], border = NA)
    }
    if (!horiz) {
      polygon(c(0, 0, 1, 1), poly[[i]], col = col[i], border = NA)
    }
  }
}