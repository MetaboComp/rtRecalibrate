#' Plot retention time (RT) adjustment
#'
#' @param raw Original data (in MsExperiment or MSnbase format)
#' @param adj RT-adjusted data (in same format as raw)
#' @param col Optional argument to specify plot color (should correspond to the order of files in raw and adj)
#' @param legend Optional argument to plot (TRUE; default) or suppress (FALSE) the legend 
#' @param cex.legend Optional argument to change font size in legend
#'
#' @return Returns a plot of RT adjustment vs the original recorded RTs
#' @export
#'
#' @examples
#' ## Read the MS data
#' mse <- MsExperiment::readMsExperiment(files)
#' 
#' ## Perform the recalibration on the data set
#' res_mse <- recalibrateRt(mse, lamasRP)
#' 
#' ## Plot RT adjustment 
#' par(mfrow = c(1, 1), mar = c(4, 4, 0, 0) + .5)
#' plotRTAdjust(mse, res_mse)
plotRTAdjust <- function(raw, adj, col, legend = TRUE, cex.legend = 0.5) {
  # Make sure relevant libraries are loaded
  if(class(raw) == 'MsExperiment') library(MsExperiment)
  if(class(raw) == 'MSnbase') library(MSnbase)
  
  # Fix plotting size for legend
  cex <- cex.legend * par()$cex
  
  # Extract filenames
  names <- basename(fileNames(raw))
  nFiles <- length(names)
  
  # Fix plotting colors
  if(missing(col)) col <- rainbow(nFiles)
  
  # Some basic sanity checks
  if(length(col) != nFiles) stop('color argument not matching length of files')
  if(any(names != basename(fileNames(adj)))) stop('File mismatch between raw and adjusted')
  
  # Extract raw and adjusted RTs
  rt_raw <- split(rtime(raw), fromFile(raw))
  rt_adj <- split(rtime(adj), fromFile(adj))
  
  # Calculate dRT from adjustments
  dRT <- rt_adj
  for (i in 1:nFiles) dRT[[i]] <- dRT[[i]] - rt_raw[[i]]
  
  # Identify plot ranges
  dRTmin <- 0
  dRTmax <- 0
  xmin <- 0
  xmax <- 0
  for (i in 1:length(dRT)) {
    dRTmin <- min(dRTmin, min(dRT[[i]]))
    dRTmax <- max(dRTmax, max(dRT[[i]]))
    xmin <- min(xmin, min(rt_raw[[i]]))
    xmax <- max(xmax, max(rt_raw[[i]]))
  }
  
  # Actual RT correction plots
  plot(0, 0, ylim = c(dRTmin, dRTmax), xlim = c(xmin, xmax), type = 'n', las = 1, 
       xlab = expression(RT[raw]), ylab = expression(RT[adj] - RT[raw]))
  for (i in 1:nFiles) {
    lines(rt_raw[[i]], 
          dRT[[i]],
          col = col[i])
  }
  
  legend('topleft', legend = names, col = col, lty = 1, bty = 'n', cex = cex)
}