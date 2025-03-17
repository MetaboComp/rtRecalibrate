#' Recalibrate retention time axis of mzML file
#'
#' @param files Path to mzML file(s)
#' @param lamas Data frame with 2 columns describing the landmarks: m/z in the first and RT in the second
#' @param method Type of warping performed: GAM or loess
#' @param dRT_roi Determines RT width of ROIs for detecting landmarks
#' @param ppm_roi Determines m/z width (in ppm) of ROIs for detecting landmarks
#' @param peakwidth CentWave parameter to XCMS for detecting peaks in ROIs
#' @param dRT_match Tolerance (dRT) for matching possible ROI peaks to landmarks
#' @param ppm_match Tolerance (dRT) for matching possible ROI peaks to landmarks
#' @param bs GAM smoothing method (defaults to thin plate; NB: B- and P-splines have been shown to produce artefacts)
#' @param span Loess smoothing parameter
#' @param ssqRatio Parameter for removing outliers in the warping
#' @param zeroWeight Adding additional weight to the zero point {0,0}, reduces "the"weird" warping in the early RT range 
#' @param jpg Whether to store diagnostic plot in jpg format in "rtRecal_log/*.jpg" (TRUE/FALSE)
#' @param plot Whether to produce plots (TRUE/FALSE; if jpg is TRUE, plots go to file instead of plot window)
#' @param save Whether to save recalibrated file in "rtRecal/*.mzML"
#' @param ... Additional arguments
#'
#' @return Depends of save, plot and jpg
#' @export
#'
#' @examples
#' rtRecal(file = file, lamas = lamasRP)
rtRecalibrate <- function(files, 
                          lamas, 
                          method = c('GAM', 'loess'),
                          dRT_roi = 15, 
                          ppm_roi = 40, 
                          dRT_match = 15, 
                          ppm_match = 5, 
                          peakwidth = c(7, 50), 
                          bs = 'tp', 
                          span = 0.5,
                          ssqRatio = 3,
                          zeroWeight = 10,
                          jpg = TRUE,
                          plot = TRUE,
                          save = TRUE,
                          ...) {
  
  library(MSnbase)
  library(xcms)
  
  corrections <- list()
  
  for (i in 1:length(files)) {
    corrections[[basename(files[i])]] <- rtRecal(file = files[i], 
                              lamas = lamas, 
                              method = method,
                              dRT_roi = dRT_roi, 
                              ppm_roi = ppm_roi, 
                              dRT_match = dRT_match, 
                              ppm_match = ppm_match, 
                              peakwidth = peakwidth, 
                              bs = bs, 
                              span = span,
                              ssqRatio = ssqRatio,
                              zeroWeight = zeroWeight,
                              jpg = jpg,
                              plot = plot,
                              save = save,
                              ...)
  }
  
  if(plot) {
    
    # Identify plot ranges
    dRTmin <- 0
    dRTmax <- 0
    xmin <- 0
    xmax <- 0
    for (i in 1:length(corrections)) {
      dRTmin <- min(dRTmin, min(corrections[[i]]$dRT))
      dRTmax <- max(dRTmax, max(corrections[[i]]$dRT))
      xmin <- min(xmin, min(corrections[[i]]$RTOld))
      xmax <- max(xmax, max(corrections[[i]]$RTOld))
    }
    
    # Colors correspondding to files
    col <- rainbow(length(corrections))
    
    # Actual RT correction plots
    plot(0, 0, ylim = c(dRTmin, dRTmax), xlim = c(xmin, xmax), type = 'n', las = 1, xlab = 'Recorded Retention Time', ylab = 'Retention Time Alignment')
    for (i in 1:length(corrections)) {
      lines(corrections[[i]]$RTOld, 
            corrections[[i]]$dRT,
            col = col[i])
    }
    
    legend('topleft', legend = names(corrections), col = col, lty = 1, bty = 'n')
  }
  
  return(invisible(corrections))
}
