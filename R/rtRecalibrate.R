#' @title Recalibrate retention times based on landmark compounds
#' @author Carl Brunius, Johannes Rainer, Anton Ribbenstedt
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
#' @param zeroWeight `numeric(1)`: additional weight to the zero point \{0,0\},
#'     reduces warping problems in the early RT range.
#'
#' @param BPPARAM Parallel processing setup. See [bpparam()] for more
#'     information.
#'     
#' @param plot Decide whether output from rtRecalibrate should be plotted or not
#'
#' @param ... ignored.
#'
#' @return Input object with recalibrated retention times.
#'
#' @import methods
#' @export
setGeneric("rtRecalibrate", def = function(object, ...)
    standardGeneric("rtRecalibrate"))

#' @importMethodsFrom Spectra rtime<-
#'
#' @importFrom BiocParallel bpparam bplapply
#' @import methods
#'
#' @examples
#' data("landmarks_final")
#' \donttest{
#'   if (requireNamespace("MsBackendMetaboLights", quietly = TRUE)) {
#'     MsBackendMetaboLights::mtbls_sync_data_files(
#'       mtblsId = "MTBLS8735",
#'       assayName = "a_MTBLS8735_LC-MS_positive_hilic_metabolite_profiling.txt"
#'     )
#'     
#'     mzMLFiles <- MsBackendMetaboLights::mtbls_cached_data_files()$rpath
#'     mzMLFiles_QC <- mzMLFiles[grepl("QC", mzMLFiles)]
#'     
#'     MsExpObj <- MsExperiment::readMsExperiment(mzMLFiles_QC)
#'     
#'     res_mse <- rtRecalibrate(object=MsExpObj,
#'                                    lamas=landmarks_final,
#'                                    method="GAM")
#'   }
#' }
#'
#' @rdname rtRecalibrate
setMethod(
    "rtRecalibrate", "MsExperiment",
    function(object, lamas, method = c("GAM", "loess"), roiExpandRt = 15,
             roiPpm = 40, toleranceRt = 15, ppm = 5, tolerance = 0,
             peakwidth = c(7, 50), bs = "tp", span = 0.5, ssqRatio = 3,
             zeroWeight = 10, BPPARAM = bpparam(), plot=TRUE) {
      
        #Error checking
        if(!(method %in% c("GAM", "loess"))){
          stop("Method has to be either \"GAM\" or \"loess\"")
        }
        if(length(method) > 1){
          stop("A single method has to be selected to use rtRecalibrate")
        }
        method <- match.arg(method)
        adj_rt <- bplapply(
            split(object, seq_along(object)),
            function(z, method, roiExpandRt, roiPpm, peakwidth, toleranceRt,
                     ppm, tolerance, bs, span, ssqRatio, zeroWeight) {
                pks <- lama_peak_detection(z, lamas, ppm = roiPpm,
                                            expandRt = roiExpandRt,
                                            peakwidth = peakwidth)
                lama_rt <- match_peaks_lamas(
                    pks, lamas, tolerance = tolerance,
                    ppm = ppm, toleranceRt = toleranceRt)
                rtm <- data.frame(reference = lamas$RT,
                                  observed = lama_rt)
                lama_rt_correction(
                    rt = rtime(z), rtMap = rtm, method = method,
                    bs = bs, span = span, ssqRatio = ssqRatio,
                    zeroWeight = zeroWeight, plot=plot)
            }, method = method, roiExpandRt = roiExpandRt, roiPpm = roiPpm,
            peakwidth = peakwidth, toleranceRt = toleranceRt, ppm = ppm,
            tolerance = tolerance, bs = bs, span = span, ssqRatio = ssqRatio,
            zeroWeight = zeroWeight, BPPARAM = BPPARAM)
        rtime(object@spectra) <- unlist(adj_rt, use.names = FALSE)
        object
    })

#' @importMethodsFrom MSnbase splitByFile fileNames
#'
#' @import methods
#' @examples
#' data("landmarks_final")
#' \donttest{
#'   if (requireNamespace("MsBackendMetaboLights", quietly = TRUE)) {
#'     MsBackendMetaboLights::mtbls_sync_data_files(
#'       mtblsId = "MTBLS8735",
#'       assayName = "a_MTBLS8735_LC-MS_positive_hilic_metabolite_profiling.txt"
#'     )
#'     
#'     mzMLFiles <- MsBackendMetaboLights::mtbls_cached_data_files()$rpath
#'     mzMLFiles_QC <- mzMLFiles[grepl("QC", mzMLFiles)]
#'     
#'     MSnExpObj <- MSnbase::readMSData(mzMLFiles_QC,
#'                                       mode="onDisk")
#'     
#'     res_mse <- rtRecalibrate(object=MSnExpObj,
#'                                    lamas=landmarks_final,
#'                                    method="GAM")
#'   }
#' }
#'
#' @rdname rtRecalibrate
setMethod(
    "rtRecalibrate", "OnDiskMSnExp",
    function(object, lamas, method = c("GAM", "loess"), roiExpandRt = 15,
             roiPpm = 40, toleranceRt = 15, ppm = 5, tolerance = 0,
             peakwidth = c(7, 50), bs = "tp", span = 0.5, ssqRatio = 3,
             zeroWeight = 10, BPPARAM = bpparam(), plot=TRUE) {
        method <- match.arg(method)
        adj_rt <- bplapply(
            splitByFile(object, f = factor(seq_along(fileNames(object)))),
            function(z, method, roiExpandRt, roiPpm, peakwidth, toleranceRt,
                     expandRt = roiExpandRt,
                     ppm, tolerance, bs, span, ssqRatio, zeroWeight) {
                pks <- lama_peak_detection(z, lamas, ppm = roiPpm,
                                            expandRt = roiExpandRt,
                                            peakwidth = peakwidth)
                lama_rt <- match_peaks_lamas(
                    pks, lamas, tolerance = tolerance,
                    ppm = ppm, toleranceRt = toleranceRt)
                rtm <- data.frame(reference = lamas$RT,
                                  observed = lama_rt)
                lama_rt_correction(
                    rt = rtime(z), rtMap = rtm, method = method,
                    bs = bs, span = span, ssqRatio = ssqRatio,
                    zeroWeight = zeroWeight, plot=plot)
            }, method = method, roiExpandRt = roiExpandRt, roiPpm = roiPpm,
            peakwidth = peakwidth, toleranceRt = toleranceRt, ppm = ppm,
            tolerance = tolerance, bs = bs, span = span, ssqRatio = ssqRatio,
            zeroWeight = zeroWeight, BPPARAM = BPPARAM)
        object@featureData$retentionTime <- unlist(adj_rt, use.names = FALSE)
        object
    })



## Add an adjustRtime,LamaParama that does the adjustment based on detected
## peaks and the lamas.
