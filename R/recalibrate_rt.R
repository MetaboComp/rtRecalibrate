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
.match_peaks_lamas <- function(chromPeaks, lamas, tolerance = 0,
                               ppm = 5, toleranceRt = 15) {
    cmz <- chromPeaks[, "mz"]
    crt <- chromPeaks[, "rt"]
    rt <- rep(NA_real_, nrow(lamas))
    for (i in seq_along(rt)) {
        res <- abs(cmz - lamas$MZ[i]) < (ppm(lamas$MZ[i], ppm) + tolerance) &
            abs(crt - lamas$RT[i]) < toleranceRt
        if (sum(res) == 1L) rt[i] <- crt[res]
    }
    rt
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
.lama_peak_detection <- function(x, lamas, ppm = 40, expandRt = 0,
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
    chromPeaks(res)
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
#'
#' @noRd
.lama_rt_correction <- function(rt, rtMap, method = c("GAM", "loess"),
                                bs = "tp", span = 0.5, ssqRatio = 3,
                                zeroWeight = 10) {
    rtMap <- rbind(rtMap, c(0, 0))
    rtMap <- rtMap[complete.cases(rtMap), ]
    rtMap <- rtMap[order(rtMap$observed), ]
    weights <- rep(1, nrow(rtMap))
    weights[1L] <- zeroWeight
    if (method == "GAM") {
        fit <- .makeGAM(RTDF = rtMap, bs = bs, ssqRatio = ssqRatio,
                        weights = weights)
    } else {
        fit <- .makeLoess(RTDF = rtMap, span = span, ssqRatio = ssqRatio,
                          weights = weights)
    }

    ## Sort out outliers and re-perform matching
    ## but make sure origo (index 1) is kept!
    fit$keep[1] <- TRUE
    rtMap <- rtMap[fit$keep, , drop = FALSE]
    weights <- weights[fit$keep]
    if (method == 'GAM') {
        fit <- .makeGAM(RTDF = rtMap, bs = bs, ssqRatio = ssqRatio,
                        weights = weights)
    } else {
        fit <- .makeLoess(RTDF = rtMap, span = span, ssqRatio = ssqRatio,
                          weights = weights)
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
    rt_new
}

#' @importFrom mgcv gam
.makeGAM <- function(RTDF, bs, ssqRatio, weights = weights) {
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

.makeLoess <- function(RTDF, span, ssqRatio, weights = weights) {
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


#' @title Recalibrate retention times based on landmark compounds
#'
#' @description
#'
#' Recalibrate retention times of provided MS data based on landmark features
#' (lamas). For input data `object` being an `MsExperiment` (from the
#' *MsExperiment* package) or `OnDiskMSnExp` (from the *MSnbase* package)
#' chromatographic peak detection will be performed specifically looking for
#' the landmark features and the retention times of the full data set will
#' be recalibrated by aligning the observed retention times of the landmark
#' features to the provided *expected* ones. Chromatographic peak detection
#' will be performed using the *centWave* method from the *xcms* package.
#'
#' Chromatographic peaks are matched to lamas based on their m/z and retention
#' times and parameters `toleranceRt`, `ppm` and `tolerance` defining
#' acceptable deviations. Only lamas/chromatographic peaks with a 1:1 mapping
#' are used to estimate retention time shifts.
#'
#' @param object Either a `MsExperiment` or `OnDiskMSnExp` object with the MS
#'     data.
#'
#' @param lamas `data.frame` with m/z values (first column) and retention
#'     times (second columns) of the landmarks.
#'
#' @param method `character(1)` with the type of warping. Either
#'     `method = "GAM"` (default) or `method = "loess"`.
#'
#' @param roiExpandRt `numeric(1)` defining the retention time width left and
#'     right of the landmark's expected retention time (second column in
#'     `lamas`) where peak detection should be performed.
#'
#' @param roiPpm `numeric(1)` defining the m/z width around the landmark's
#'     expected m/z (first column in `lamas`) where peak detection should be
#'     performed.
#'
#' @param toleranceRt `numeric(1)` with the maximal allowed (absolute)
#'     difference in retention times between lamas and chromatographic peaks.
#'     Used for the mapping of identified chromatographic peaks and lamas.
#'
#' @param ppm `numeric(1)` defining the m/z-relative maximal allowed
#'     difference in m/z between lamas and chromatographic peaks. Used for
#'     the mapping of identified chromatographic peaks and lamas.
#'
#' @param tolerance `numeric(1)` defining the absolute acceptable difference
#'     in m/z between lamas and chromatographic peaks. Used for the mapping
#'     of identified chromatographic peaks and lamas.
#'
#' @param peakwidth `numeric(2)` with the lower and upper limit (in retention
#'     time) for the *centWave*-based chromatographic peak detection.
#'
#' @param bs `character(1)` defining the GAM moothing method. (defaults to
#'     thin plate; NB: B- and P-splines have been shown to produce artefacts).
#'
#' @param span `numeric(1)` loess smoothing parameter.
#'
#' @param ssqRatio `numeric(1)` defining the settings for removing outliers
#'     in the warping.
#'
#' @param zeroWeight `numeric(1)`: additional weight to the zero point {0,0},
#'     reduces warping problems in the early RT range.
#'
#' @param BPPARAM Parallel processing setup. See [bpparam()] for more
#'     information.
#'
#' @param ... ignored.
#'
#' @return Input object with recalibrated retention times.
#'
#' @author Carl Brunius, Johannes Rainer.
#'
#' @exportMethod recalibrateRt
setGeneric("recalibrateRt", def = function(object, ...)
    standardGeneric("recalibrateRt"))

#' @importMethodsFrom Spectra rtime<-
#'
#' @importFrom BiocParallel bpparam bplapply
#'
#' @rdname recalibrateRt
setMethod(
    "recalibrateRt", "MsExperiment",
    function(object, lamas, method = c("GAM", "loess"), roiExpandRt = 15,
             roiPpm = 40, toleranceRt = 15, ppm = 5, tolerance = 0,
             peakwidth = c(7, 50), bs = "tp", span = 0.5, ssqRatio = 3,
             zeroWeight = 10, BPPARAM = bpparam()) {
        method <- match.arg(method)
        adj_rt <- bplapply(
            split(object, seq_along(object)),
            function(z, method, roiExpandRt, roiPpm, peakwidth, toleranceRt,
                     ppm, tolerance, bs, span, ssqRatio, zeroWeight) {
                pks <- .lama_peak_detection(z, lamas, ppm = roiPpm,
                                            expandRt = roiExpandRt,
                                            peakwidth = peakwidth)
                lama_rt <- .match_peaks_lamas(
                    pks, lamas, tolerance = tolerance,
                    ppm = ppm, toleranceRt = toleranceRt)
                rtm <- data.frame(reference = lamas$RT,
                                  observed = lama_rt)
                .lama_rt_correction(
                    rt = rtime(z), rtMap = rtm, method = method,
                    bs = bs, span = span, ssqRatio = ssqRatio,
                    zeroWeight = zeroWeight)
            }, method = method, roiExpandRt = roiExpandRt, roiPpm = roiPpm,
            peakwidth = peakwidth, toleranceRt = toleranceRt, ppm = ppm,
            tolerance = tolerance, bs = bs, span = span, ssqRatio = ssqRatio,
            zeroWeight = zeroWeight, BPPARAM = BPPARAM)
        rtime(object@spectra) <- unlist(adj_rt, use.names = FALSE)
        object
    })

#' @importMethodsFrom MSnbase splitByFile fileNames
#'
#' @rdname recalibrateRt
setMethod(
    "recalibrateRt", "OnDiskMSnExp",
    function(object, lamas, method = c("GAM", "loess"), roiExpandRt = 15,
             roiPpm = 40, toleranceRt = 15, ppm = 5, tolerance = 0,
             peakwidth = c(7, 50), bs = "tp", span = 0.5, ssqRatio = 3,
             zeroWeight = 10, BPPARAM = bpparam()) {
        method <- match.arg(method)
        adj_rt <- bplapply(
            splitByFile(object, f = factor(seq_along(fileNames(object)))),
            function(z, method, roiExpandRt, roiPpm, peakwidth, toleranceRt,
                     ppm, tolerance, bs, span, ssqRatio, zeroWeight) {
                pks <- .lama_peak_detection(z, lamas, ppm = roiPpm,
                                            expandRt = roiExpandRt,
                                            peakwidth = peakwidth)
                lama_rt <- .match_peaks_lamas(
                    pks, lamas, tolerance = tolerance,
                    ppm = ppm, toleranceRt = toleranceRt)
                rtm <- data.frame(reference = lamas$RT,
                                  observed = lama_rt)
                .lama_rt_correction(
                    rt = rtime(z), rtMap = rtm, method = method,
                    bs = bs, span = span, ssqRatio = ssqRatio,
                    zeroWeight = zeroWeight)
            }, method = method, roiExpandRt = roiExpandRt, roiPpm = roiPpm,
            peakwidth = peakwidth, toleranceRt = toleranceRt, ppm = ppm,
            tolerance = tolerance, bs = bs, span = span, ssqRatio = ssqRatio,
            zeroWeight = zeroWeight, BPPARAM = BPPARAM)
        object@featureData$retentionTime <- unlist(adj_rt, use.names = FALSE)
        object
    })



## Add an adjustRtime,LamaParama that does the adjustment based on detected
## peaks and the lamas.
