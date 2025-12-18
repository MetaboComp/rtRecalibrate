#' Recalibrate retention time axis of mzML file to landmark reference
#'  - One file only
#'
#' Please see rtRecalibrate for standard interface for
#'  multiple files and parallelization
#'
#' @param file Path to a mzML file
#' @param lamas Data frame with 2 columns describing the landmarks: \"MZ\" in
#'     the first and \"RT\" in the second
#' @param method Type of warping performed: GAM or loess
#' @param dRT_roi Determines RT width of ROIs for detecting landmarks
#' @param ppm_roi Determines m/z width (in ppm) of ROIs for detecting landmarks
#' @param dRT_match Tolerance (dRT) for matching possible ROI peaks to landmarks
#' @param ppm_match Tolerance (dRT) for matching possible ROI peaks to landmarks
#' @param peakwidth CentWave parameter to XCMS for detecting peaks in ROIs
#' @param expandRt Defines the retention time range of the
#'     ROIs in lama_peak_detection(). The retention time range will be the
#'     expected retention time of a lama +/- `expandRt`.
#' @param bs GAM smoothing method (defaults to thin plate; NB: B- and P-splines
#'      have been shown to produce artefacts)
#' @param span Loess smoothing parameter
#' @param ssqRatio Parameter for removing outliers in the warping
#' @param zeroWeight Adding additional weight to the zero point \{0,0\},
#'      reduces "the"weird" warping in the early RT range
#' @param jpg Whether to store diagnostic plot in jpg format in
#'      "rtRecal_log/*.jpg" (TRUE/FALSE)
#' @param plot Whether to produce plots (TRUE/FALSE; if jpg is TRUE, plots go
#'      to file instead of plot window)
#' @param save Whether to save recalibrated file in "rtRecal/*.mzML"
#'
#' @return Depends of save, plot and jpg
#' @importFrom grDevices jpeg
#' @importFrom graphics abline
#' @importFrom graphics legend
#' @importFrom graphics hist
#' @importFrom grDevices dev.off
#' @importFrom xcms writeMSData
#'

rtRecal <- function(file,
                    lamas,
                    method = c("GAM", "loess"),
                    dRT_roi = 15,
                    ppm_roi = 40,
                    dRT_match = 15,
                    ppm_match = 5,
                    peakwidth = c(7, 50),
                    expandRt = 0,
                    bs = "tp",
                    span = 0.5,
                    ssqRatio = 3,
                    zeroWeight = 10,
                    jpg = TRUE,
                    plot = TRUE,
                    save = TRUE) {
    ## Testing if file exists (implicitly checks directory as well)
    if (!file.exists(file)) {
        stop("The file doesn't exist. Please choose the path to a real file.")
    }

    ## Testing that file is mzML format
    if (!grepl(".mzML", file)) {
        stop("The file chosen ist not a .mzML file.")
    }

    ## Testing that user has chosen one of the two methods
    if (length(method) > 1) {
        stop("Choose one of the two methods for the alignment.")
    }

    if (nrow(lamas) < 5 || ncol(lamas) > 2 || !is.data.frame(lamas)) {
        stop("Wrong format for lamas input. Should be a data frame with
        2 columns \"mz\" and \"rt\" and with at least 5 lamas")
    }

    ## Setup storage
    if (jpg) {
        plot <- TRUE
        if (!dir.exists(paste0(dirname(file), "/rtRecal_log"))) {
            dir.create(paste0(dirname(file), "/rtRecal_log"))
        }
        if (file.exists(paste0(
            dirname(file), "/rtRecal_log/",
            basename(file), ".jpg"
        ))) {
            file.remove(paste0(dirname(file),
                        "/rtRecal_log/", basename(file), ".jpg"))
        }

        jpeg(
            filename = paste0(
                dirname(file), "/rtRecal_log/",
                basename(file), ".jpg"
            ), width = 1000,
            height = 2000, pointsize = 30
        )
    }

    if (save) {
        if (!dir.exists(paste0(dirname(file), "/rtRecal"))) {
            dir.create(paste0(dirname(file), "/rtRecal"))
        }
        if (file.exists(paste0(dirname(file), "/rtRecal/", basename(file)))) {
            file.remove(paste0(dirname(file), "/rtRecal/", basename(file)))
        }
    }

    ## Reading raw data and detecting landmarks in files
    ms <- MSnbase::readMSData(file,
        msLevel = 1,
        verbose = FALSE,
        mode = "onDisk"
    )
    rt <- MSnbase::rtime(ms)

    peaks <- lama_peak_detection(
        x = ms,
        lamas = lamas,
        ppm = ppm_roi,
        expandRt = expandRt,
        peakwidth = peakwidth
    )
    peaks <- as.data.frame(peaks)

    ## Restrict found peaks to (more-or-less) unambiguous lamas
    matchRT <- match_peaks_lamas(peaks,
        lamas,
        ppm = ppm_match,
        toleranceRt = dRT_match
    )

    ## Combine observed (sample) and reference (lama) times in a dataframe
    rtMap <- data.frame(reference = lamas$RT, observed = matchRT)

    ## Correct retention times in files based on reference RTs
    rt_new <- lama_rt_correction(rt,
        rtMap = rtMap,
        method = method,
        bs = bs,
        span = span,
        ssqRatio = ssqRatio,
        zeroWeight = zeroWeight,
        plot = plot
    )

    ## Turn off development of jpgs
    if (jpg) dev.off()

    ## Save file if specified by user
    if (save) {
        ms@featureData$retentionTime <- rt_new
        writeMSData(ms, file = paste0(
            dirname(file),
            "/rtRecal/", basename(file)
        ), copy = TRUE)
    }

    return(invisible(data.frame(
        RTOld = rt,
        RTNew = rt_new,
        dRT = rt_new - rt
    )))
}
