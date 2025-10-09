#' @title Match identified chromatographic peaks to lamas
#'
#' @description
#'
#' Given a matrix with identified chromatographic peaks, identify those that
#' match (given RT and m/z tolerances) to those of lamas and report their
#' retention times.
#'
#' @param chromPeaks `matrix` with identified chromatographic peaks of **one**
#'     sample.
#'
#' @param lamas `data.frame` with retention times (`$RT`) and m/z (`$MZ`) of
#'     lamas.
#'
#' @param tolerance `numeric(1)` with the absolute difference in m/z values.
#'
#' @param ppm `numeric(1)` with the m/z-relative tolerance.
#'
#' @param toleranceRt `numeric(1)` with the absolute acceptable difference in
#'     retention times.
#'
#' @return `numeric` with retention times.
#'
#' @noRd
#'
#' @importFrom MsCoreUtils ppm
match_peaks_lamas <- function(chromPeaks,
                              lamas, tolerance = 0,
                               ppm = 5, toleranceRt = 15) {
  cmz <- chromPeaks[, "mz"]
  crt <- chromPeaks[, "rt"]
  rt <- rep(NA_real_, nrow(lamas))
  for (i in seq_along(rt)) {
    res <- abs(cmz - lamas$MZ[i]) < (ppm(lamas$MZ[i], ppm) + tolerance) &
      abs(crt - lamas$RT[i]) < toleranceRt
    if (sum(res) == 1L) rt[i] <- crt[res]
  }
  return(rt)
}

#' @title Perform peak detection for lamas in a single data file
#'
#' @param x `OnDiskMSnExp` (*MSnbase*) or `MsExperiment` (*MsExperiment*) object
#'     (single file)
#'
#' @param lamas `data.frame` with the LaMas
#'
#' @param ppm `numeric(1)` defining the ppm by which the m/z range of the
#'     regions of interest (ROI) for each lama should be defined.
#'
#' @param expandRt `numeric(1)` used to define the retention time range of the
#'     ROIs. The retention time range will be the expected retention time
#'     of a lama +/- `expandRt`.
#'
#' @param peakwidth `numeric(2)` defining the lower and upper peak width (in
#'     retention time dimension) used for *centWave*-based peak detection.
#'
#' @return `matrix` with identified chromatographic peaks.
#'
#' @importMethodsFrom xcms rtime findChromPeaks chromPeaks
#'
#' @importFrom xcms CentWaveParam
#'
#' @importFrom BiocParallel SerialParam
#'
#' @noRd
lama_peak_detection <- function(x, lamas, ppm = 40, expandRt = 0,
                                 peakwidth = c(7, 50)) {
  rt_raw <- rtime(x)
  ## Define ROIs.
  rl <- vector("list", nrow(lamas))
  for (i in seq_len(nrow(lamas))) {
    roi <- list()
    roi$mz <- lamas$MZ[i]
    roi$mzmin <- roi$mz - ppm(roi$mz, ppm)
    roi$mzmax <- roi$mz + ppm(roi$mz, ppm)
    roi$scmin <- which.min(abs(rt_raw - (lamas$RT[i] - expandRt)))
    roi$scmax <- which.min(abs(rt_raw - (lamas$RT[i] + expandRt)))
    roi$length <- roi$scmax - (roi$scmin - 1)
    roi$intensity <- 1e5 # maybe provide that as a parameter?
    rl[[i]] <- roi
  }
  ## Perform peak detection.
  cwp <- CentWaveParam(peakwidth = peakwidth, roiList = rl)
  res <- findChromPeaks(x, param = cwp, BPPARAM = SerialParam())
  return(chromPeaks(res))
}

#' @title Base LaMa-based retention correction function
#'
#' @description
#'
#' Adjusts provided retention times **of a single data file** using observed
#' and theoretical retention times of LaMas.
#'
#' @param rt `numeric` with the raw retention times of a sample.
#'
#' @param rtMap `data.frame` with *reference* retention times of LaMas and
#'     *observed* retention times of matching peaks in the same sample from
#'     which the retention times in `rt` are.
#'
#' @return `numeric` with the adjusted retention times
#'
#' @importFrom stats complete.cases predict resid
#' @importFrom grDevices dev.new
#'
#' @noRd
lama_rt_correction <- function(rt,
                               rtMap,
                               method = c("GAM", "loess"),
                               bs = "tp",
                               span = 0.5,
                               ssqRatio = 3,
                               zeroWeight = 10,
                               plot) {
  rtMap <- rbind(rtMap, c(0, 0))
  rtMap <- rtMap[complete.cases(rtMap), ]
  rtMap <- rtMap[order(rtMap$observed), ]
  weights <- rep(1, nrow(rtMap))
  weights[1L] <- zeroWeight
  if (method == "GAM") {
    fit <- makeGAM(RTDF = rtMap, bs = bs, ssqRatio = ssqRatio,
                    weights = weights)
  } else {
    fit <- makeLoess(RTDF = rtMap, span = span, ssqRatio = ssqRatio,
                      weights = weights)
  }
  
  ## Sort out outliers and re-perform matching
  ## but make sure origo (index 1) is kept!
  fit$keep[1] <- TRUE
  rtMap <- rtMap[fit$keep, , drop = FALSE]
  weights <- weights[fit$keep]
  
  ## Preparing plot
  if(plot){
    dev.new()
    par(mar = c(4, 4, 0, 0) + .5, fig = c(0, 1, 0.5, 1))
  } 
  
  if (method == 'GAM') {
    fit <- makeGAM(RTDF = rtMap, bs = bs, ssqRatio = ssqRatio,
                    weights = weights)
  } else {
    fit <- makeLoess(RTDF = rtMap, span = span, ssqRatio = ssqRatio,
                      weights = weights)
  }
  ## Making plot of fit 
  if(plot) {
    plot((reference - observed) ~ observed, data = rtMap,
         col = ifelse(fit$keep, 'black', 'red'), xlim = c(0, max(rt)),
         xlab = expression(rt[raw]), ylab = expression(rt[adj] - rt[raw]))
    abline(h = 0, col = 'grey', lty = 2)
    # lines(rtDF$observed, fit$pred - rtDF$observed)
  }
  
  
  ## Apply last corrected scan (i.e. at last observed lama) to all
  ## subsequent scans
  dRTEnd <- (fit$pred - rtMap$observed)[nrow(rtMap)]
  rt_new <- rt
  rt_mask <- rt <= rtMap$observed[nrow(rtMap)]
  rt_new[rt_mask] <- predict(fit$model,
                             newdata = data.frame(observed = rt[rt_mask]))
  rt_mask <- rt > rtMap$observed[nrow(rtMap)]
  rt_new[rt_mask] <- rt[rt_mask] + dRTEnd
  rt_new[rt_new < 0] <- 0
  
  ## Plot RT correction applied to all scans
  if(plot) {
    lines(rt, rt_new - rt, type = 'l',
          main = 'Final correction all scans - RT-diff')
    
    legend('topleft', pch = c(1, 1, NA), lty = c(NA, NA, 1),
           col = c('black', 'red', 'black'),
           legend = c('LaMa (used)', 'LaMa (excluded)', 'RT adjustment'),
           bty = 'n', cex = 0.75 * par()$cex)
  }
  
  
  
  # Plot histograms of between-scan dRT
  if(plot) {
    dRT_org <- rt[-1] - rt[-length(rt)]
    dRT_corr <- rt_new[-1] - rt_new[-length(rt_new)]
    par(fig = c(0, 0.5, 0, 0.5), new = TRUE)
    hist(dRT_org, xlim = range(c(dRT_org, dRT_corr)), 
         xlab = expression(dRT[Raw]), main = '')
    par(fig = c(0.5, 1, 0, 0.5), new = TRUE)
    hist(dRT_corr, xlim = range(c(dRT_org, dRT_corr)), 
         xlab = expression(dRT[Adj]), main = '')
  }
  
  return(rt_new)
}

#' @importFrom mgcv gam
makeGAM <- function(RTDF, bs, ssqRatio, weights = weights) {
  gam <- gam(reference ~ s(observed, bs = bs), weights = weights,
             data = RTDF)
  pred <- predict(gam)
  SSq <- resid(gam)^2
  meanSSq <- mean(SSq)
  keep <- ifelse((SSq / meanSSq) > ssqRatio, FALSE, TRUE)
  return(list(model = gam,
              pred = pred,
              keep = keep))
}

makeLoess <- function(RTDF, span, ssqRatio, weights = weights) {
  loess <- loess(reference ~ observed, weights = weights, span = span,
                 data = RTDF)
  pred <- predict(loess)
  SSq <- resid(loess)^2
  meanSSq <- mean(SSq)
  keep <- ifelse((SSq / meanSSq) > ssqRatio, FALSE, TRUE)
  return(list(model = loess,
              pred = pred,
              keep = keep))
}