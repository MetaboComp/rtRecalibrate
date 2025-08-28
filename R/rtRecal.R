#' Recalibrate retention time axis of mzML file to landmark reference
#'  - One file only
#'
#' Please see rtRecalibrate for standard interface for
#'  multiple files and parallelization
#'
#' @param file Path to a mzML file
#' @param lamas Data frame with 2 columns describing the landmarks: m/z in the
#'  first and RT in the second
#' @param method Type of warping performed: GAM or loess
#' @param dRT_roi Determines RT width of ROIs for detecting landmarks
#' @param ppm_roi Determines m/z width (in ppm) of ROIs for detecting landmarks
#' @param peakwidth CentWave parameter to XCMS for detecting peaks in ROIs
#' @param dRT_match Tolerance (dRT) for matching possible ROI peaks to landmarks
#' @param ppm_match Tolerance (dRT) for matching possible ROI peaks to landmarks
#' @param bs GAM smoothing method (defaults to thin plate; NB: B- and P-splines 
#' have been shown to produce artefacts)
#' @param span Loess smoothing parameter
#' @param ssqRatio Parameter for removing outliers in the warping
#' @param zeroWeight Adding additional weight to the zero point {0,0}, 
#' reduces "the"weird" warping in the early RT range
#' @param jpg Whether to store diagnostic plot in jpg format in 
#' "rtRecal_log/*.jpg" (TRUE/FALSE)
#' @param plot Whether to produce plots (TRUE/FALSE; if jpg is TRUE, plots go 
#' to file instead of plot window)
#' @param save Whether to save recalibrated file in "rtRecal/*.mzML"
#'
#' @return Depends of save, plot and jpg
#' @export
#'
#' @examples
#' rtRecal(file = file, lamas = lamasRP)
rtRecal <- function(file,
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
                    save = TRUE) {

  # Setup storage
  if (jpg) {
    plot <- TRUE
    if(!dir.exists(paste0(dirname(file), '/rtRecal_log'))){
      dir.create(paste0(dirname(file), '/rtRecal_log'))
    }
    if(file.exists(paste0(dirname(file), '/rtRecal_log/',
                          basename(file),'.jpg'))){
      file.remove(paste0(dirname(file), '/rtRecal_log/', basename(file),'.jpg'))
    }
      
    jpeg(filename = paste0(dirname(file), '/rtRecal_log/',
                           basename(file),'.jpg'), width = 1000,
                           height = 2000, pointsize = 30)
  }

  if(save) {
    if(!dir.exists(paste0(dirname(file), '/rtRecal'))){
      dir.create(paste0(dirname(file), '/rtRecal'))
    }
    if(file.exists(paste0(dirname(file), '/rtRecal/', basename(file)))){
      file.remove(paste0(dirname(file), '/rtRecal/', basename(file)))
    }
  }

  if(missing(method)) method <- 'GAM'

  # Read raw data and extract header info
  MS <- MSnbase::readMSData(file, msLevel = 1, verbose = FALSE, mode = 'onDisk')
  # MS <- readMSData(file, msLevel = 1, verbose = FALSE)
  nScan <- length(MS)
  RT <- MSnbase::rtime(MS)

  # Convert LaMas to ROIs for xcms
  roiList <-  list()
  for (r in 1:nrow(lamas)) {
    roi <- list()
    roi$mz <- lamas[r, ]$MZ %>% unlist
    roi$mzmin <- roi$mz - ppm_roi * 1e-6 * roi$mz %>% unlist
    roi$mzmax <- roi$mz + ppm_roi * 1e-6 * roi$mz %>% unlist
    roi$scmin <- which.min(abs(RT - (lamas[r, ]$RT - dRT_roi))) %>% unlist
    roi$scmax <- which.min(abs(RT - (lamas[r, ]$RT + dRT_roi))) %>% unlist
    roi$length <- roi$scmax - (roi$scmin - 1) %>% unlist
    roi$intensity <- 1e5 %>% unlist
    roiList[[r]] <- roi
  }

  # Use XCMS to pick peaks in lama-regions
  cwParam <- xcms::CentWaveParam(roiList = roiList, peakwidth = peakwidth)
  lamaFind <- xcms::findChromPeaks(object = MS, param = cwParam)
  peaks <- xcms::chromPeaks(lamaFind) %>% as.data.frame()

  # Restrict found peaks to (more-or-less) unambiguous lamas
  # i.e. within a dppm window and only 1 peak in the RT-range
  matchRT <- rep(NA, nrow(lamas))
  matches <- rep(NA, nrow(lamas))
  for (r in 1:nrow(lamas)) {
    mz <- lamas[r, ]$MZ
    rt <- lamas[r, ]$RT
    match <- (abs(mz - peaks$mz) < (ppm_match * 1e-6 * mz)) &
                (abs(rt - peaks$rt) < dRT_match)
    matches[r] <- sum(match)
    if (sum(match) == 1) matchRT[r] <- peaks$rt[match]
  }

  # Combine observed (sample) and reference (lama) times in a dataframe
  rtDF <- data.frame(reference = lamas$RT, observed = matchRT)
  rtDF[nrow(rtDF) + 1,] <- 0
  rtDF <- rtDF[complete.cases(rtDF),]
  rtDF <- rtDF[order(rtDF$observed),]
  weights <- rep(1, nrow(rtDF))
  weights[1] <- zeroWeight

  makeGAM <- function(RTDF, bs, ssqRatio, weights = weights) {
    gam <- mgcv::gam(reference ~ s(observed, bs = bs),
                     weights = weights, data = RTDF)
    pred <- predict(gam)
    SSq <- resid(gam)^2
    meanSSq <- mean(SSq)
    keep <- ifelse((SSq / meanSSq) > ssqRatio, FALSE, TRUE)
    return(list(model = gam,
                pred = pred,
                keep = keep))
  }

  makeLoess <- function(RTDF, span, ssqRatio, weights = weights) {
    loess <- loess(reference ~ observed, weights = weights,
                   span = span, data = RTDF)
    pred <- predict(loess)
    SSq <- resid(loess)^2
    meanSSq <- mean(SSq)
    keep <- ifelse((SSq / meanSSq) > ssqRatio, FALSE, TRUE)
    return(list(model = loess,
                pred = pred,
                keep = keep))
  }

  # Prepare for plots
  if(plot) par(mar = c(4, 4, 0, 0) + .5, fig = c(0, 1, 0.5, 1))

  # par(mfrow = c(2, 3), mar = c(4, 4, 2, 0) + .5)

  # Make 1st round of RT matching
  if(method == 'GAM') {
    fit <- makeGAM(RTDF = rtDF, bs = bs, ssqRatio = ssqRatio, weights = weights)
  } else {
    fit <- makeLoess(RTDF = rtDF, span = span,
                     ssqRatio = ssqRatio, weights = weights)
  }
  if(plot) {
    plot((reference - observed) ~ observed, data = rtDF,
         col = ifelse(fit$keep, 'black', 'red'), xlim = c(0, max(RT)),
         xlab = expression(RT[raw]), ylab = expression(RT[adj] - RT[raw]))
    abline(h = 0, col = 'grey', lty = 2)
    # lines(rtDF$observed, fit$pred - rtDF$observed)
  }

  # Sort out outliers and re-perform matching
  # but make sure origo (index 1) is kept!
  fit$keep[1] <- TRUE
  rtDF <- rtDF[fit$keep,]
  weights <- weights[fit$keep]

  if(method == 'GAM') {
    fit <- makeGAM(RTDF = rtDF, bs = bs, ssqRatio = ssqRatio, weights = weights)
  } else {
    fit <- makeLoess(RTDF = rtDF, span = span,
                     ssqRatio = ssqRatio, weights = weights)
  }

  # Apply last corrected scan (i.e. at last observed lama) to subsequent scans
  dRTEnd <- (fit$pred - rtDF$observed)[nrow(rtDF)]
  RTNew <- RT
  RTNew[RT <= rtDF$observed[nrow(rtDF)]] <- 
    predict(fit$model,
            newdata = data.frame(observed =
                                   RTNew[RT <= rtDF$observed[nrow(rtDF)]]))
  RTNew[RT > rtDF$observed[nrow(rtDF)]] <- 
    RT[RT > rtDF$observed[nrow(rtDF)]] + dRTEnd
  
  #Remove all < 0
  RTNew[RTNew < 0] <- 0

  # Plot RT correction applied to all scans
  if(plot) {
    lines(RT, RTNew - RT, type = 'l',
          main = 'Final correction all scans - RT-diff')
    # abline(h = 0, col = 'grey', lty = 2)
  }
  
  legend('topleft', pch = c(1, 1, NA), lty = c(NA, NA, 1),
         col = c('black', 'red', 'black'),
         legend = c('LaMa (used)', 'LaMa (excluded)', 'RT adjustment'),
         bty = 'n', cex = 0.75 * par()$cex)

  # Plot histograms of between-scan dRT
  if(plot) {
    dRT_org <- RT[-1] - RT[-length(RT)]
    dRT_corr <- RTNew[-1] - RTNew[-length(RTNew)]
    par(fig = c(0, 0.5, 0, 0.5), new = TRUE)
    hist(dRT_org, xlim = range(c(dRT_org, dRT_corr)), 
         xlab = expression(dRT[Raw]), main = '')
    par(fig = c(0.5, 1, 0, 0.5), new = TRUE)
    hist(dRT_corr, xlim = range(c(dRT_org, dRT_corr)), 
         xlab = expression(dRT[Adj]), main = '')
  }

  if(jpg) dev.off()

  if(save) {
    MS@featureData$retentionTime <- RTNew
    writeMSData(MS, file = paste0(dirname(file),
                                  '/rtRecal/', basename(file)), copy = TRUE)
  }

  return(invisible(data.frame(RTOld = RT,
                              RTNew = RTNew,
                              dRT = RTNew - RT)))

}
