#' optimizeFindLM - A function for optimization of finding new LaMas
#'
#' @param PT_filepath File path to a .rds file which contains a peak table
#' generated from an XCMS object which has been processed until,
#' but not including, the gap-filling step
#' @param PT A peak table generated from an XCMS object which has been
#' processed until, but not including, the gap-filling step
#' @param split The character used in the XCMS object for splitting up mz & RT
#' @param minrttocheck Dead volume cut-off in beginning of chromatogram to
#' disregard when looking for LaMas
#' @param Prefilterintensity Prefilter intensity settings used when peak
#' picking with XCMS
#' @param minLM The minimum number of landmarks which are required to be
#' returned by the optimization in order to count as a result
#' @param maxLM The maxmimum number of landmarks allowed to be included in
#' the optimization
#' @return Returning a list containing 10 sets of optimized LaMas, the result
#' parameters and a df of the result parameters
#' @export
#' @examples
#' data("peak_table")
#' library(rtRecalibrate)
#'
#' landmark_suggestions <- optimizeFindLM(
#'     PT = peak_table,
#'     split = "@",
#'     Prefilterintensity = 7000,
#'     minLM = 40
#' )
#'
#' landmarks_final <- landmark_suggestions[[1]][[2]]
#'
#' print(landmarks_final)
optimizeFindLM <- function(PT_filepath = NULL,
                            PT = NULL,
                            split = "@",
                            minrttocheck = 40,
                            Prefilterintensity,
                            minLM = 100,
                            maxLM = 200) {
    ## Error checking input variables
    if (!missing(PT_filepath) && !missing(PT)) {
        stop("Cannot determine which to use if both path and object supplied!
            Choose one and remove the other.")
    }

    if (!missing(PT_filepath)) {
        if (!file.exists(PT_filepath)) {
            stop("The file doesn't exist!")
        }

        if (!any(grepl(".mzML", PT_filepath))) {
            stop("")
        }
        PT_NoFill <- readRDS(PT_filepath)
    } else if (!missing(PT)) {
        PT_NoFill <- PT
    } else {
        stop("Neither peak table nor a filepath for .rds containing peak table
        supplied!")
    }

    if (missing(Prefilterintensity)) {
        stop("The prefilter intensity (xcms) or equivalent (in other peak 
        picking software) needs to be supplied in order to correctly
        determine which peaks are eligble to become landmarks.")
    }

    if (ncol(PT_NoFill) < minLM) {
        stop("The number of features in the peak table are not enough to 
        satisfy the minLM criterium. Please lower the minLM or look over
        your peak table.")
    }

    preFCPnames <- colnames(PT_NoFill)
    splitnames <- strsplit(preFCPnames, split = split)
    mzvec <- as.numeric(vapply(splitnames, "[[", character(1), 1))
    rtvec <- as.numeric(vapply(splitnames, "[[", character(1), 2))
    runTime <- max(rtvec)
    runTime <- round(runTime, digits = 0)
    mzrtdf <- data.frame(mzvec, rtvec)
    minrttocheck <- minrttocheck
    minIntensity <- seq(
        from = Prefilterintensity,
        to = Prefilterintensity * 5,
        length.out = 10
    )
    allowedmissingness <- c(0.25, 0.20, 0.15, 0.1, 0.05, 0.01)
    mzdiff <- seq(from = 0.1, to = 1, by = 0.1)
    rtdiff <- seq(from = runTime / 36, to = runTime / 18, length.out = 8)
    Des <- expand.grid(minIntensity, allowedmissingness, mzdiff, rtdiff)
    colnames(Des) <- c("minIntensity", "allowedmissingness", "mzdiff", "rtdiff")
    foundLMs <- list()
    for (i in seq_len(nrow(Des))) {
        foundLMs[[i]] <- findLandmarks(
            dat = PT_NoFill,
            mzrtdf = mzrtdf,
            mzdiff = Des$mzdiff[i],
            rtdiff = Des$rtdiff[i],
            minIntensity = Des$minIntensity[i],
            minrttocheck = minrttocheck,
            allowedmissingness = Des$allowedmissingness[i]
        )
    }

    nLMvec <- vapply(foundLMs, FUN = nrow, FUN.VALUE = numeric(1))
    whcorrectnLM <- which(nLMvec < maxLM & nLMvec > minLM)
    Coverage <- vector()
    Settings <- Des[whcorrectnLM, ]
    for (i in seq_len(length(whcorrectnLM))) {
        Settingcandidate <- foundLMs[[whcorrectnLM[i]]]
        Coverage[i] <- calcCoverage(Settingcandidate[, 2], runTime)
    }
    Resmatr <- cbind(Settings, Coverage)
    range01 <- function(x) {
        (x - min(x)) / (max(x) - min(x))
    }
    Scaledmatr <- apply(Resmatr, 2, FUN = range01)
    if (any(is.nan(Scaledmatr))) {
        Scaledmatr[is.nan(Scaledmatr)] <- 1
    }
    Scaledmatr[, 2] <- Scaledmatr[, 2] * -1
    testobjfun <- apply(Scaledmatr, 1, sum)
    sortedscore <- sort(testobjfun, decreasing = TRUE)
    sortedscore <- as.numeric(names(sortedscore)[seq_len(10)])
    sortedCoverage <- order(testobjfun, decreasing = TRUE)
    rank1 <- foundLMs[[sortedscore[1]]]
    rank1settings <- cbind(Des[sortedscore[1], ], Coverage[sortedCoverage[1]])
    rank2 <- foundLMs[[sortedscore[2]]]
    rank2settings <- cbind(Des[sortedscore[2], ], Coverage[sortedCoverage[2]])
    rank3 <- foundLMs[[sortedscore[3]]]
    rank3settings <- cbind(Des[sortedscore[3], ], Coverage[sortedCoverage[3]])
    rank4 <- foundLMs[[sortedscore[4]]]
    rank4settings <- cbind(Des[sortedscore[4], ], Coverage[sortedCoverage[4]])
    rank5 <- foundLMs[[sortedscore[5]]]
    rank5settings <- cbind(Des[sortedscore[5], ], Coverage[sortedCoverage[5]])
    rank6 <- foundLMs[[sortedscore[6]]]
    rank6settings <- cbind(Des[sortedscore[6], ], Coverage[sortedCoverage[6]])
    rank7 <- foundLMs[[sortedscore[7]]]
    rank7settings <- cbind(Des[sortedscore[7], ], Coverage[sortedCoverage[7]])
    rank8 <- foundLMs[[sortedscore[8]]]
    rank8settings <- cbind(Des[sortedscore[8], ], Coverage[sortedCoverage[8]])
    rank9 <- foundLMs[[sortedscore[9]]]
    rank9settings <- cbind(Des[sortedscore[9], ], Coverage[sortedCoverage[9]])
    rank10 <- foundLMs[[sortedscore[10]]]
    rank10settings <- cbind(Des[sortedscore[10], ], 
                            Coverage[sortedCoverage[10]])
    LMs <- list(
        rank1,
        rank2,
        rank3,
        rank4,
        rank5,
        rank6,
        rank7,
        rank8,
        rank9,
        rank10
    )
    Settings <- list(
        rank1settings,
        rank2settings,
        rank3settings,
        rank4settings,
        rank5settings,
        rank6settings,
        rank7settings,
        rank8settings,
        rank9settings,
        rank10settings
    )
    Overview <- rbind(
        as.matrix(rank1settings),
        as.matrix(rank2settings),
        as.matrix(rank3settings),
        as.matrix(rank4settings),
        as.matrix(rank5settings),
        as.matrix(rank6settings),
        as.matrix(rank7settings),
        as.matrix(rank8settings),
        as.matrix(rank9settings),
        as.matrix(rank10settings)
    )
    colnames(Overview) <- c(
        "minIntensity",
        "allowedMissingness",
        "mzdiff",
        "rtdiff",
        "coverage"
    )
    Out <- list(LMs, Settings, Overview)
    return(Out)
}

makeInterval <- function(rts) {
    rts <- round(rts, digits = 0)
    interval <- seq(from = rts - 20, to = rts + 20, by = 1)

    return(interval)
}

calcCoverage <- function(rts, runTime) {
    interval <- unlist(lapply(rts, makeInterval), use.names = FALSE)

    interval <- interval[which(!interval > runTime)]
    interval <- interval[which(!interval < 0)]
    coverage <- length(unique(interval)) / (runTime + 1)

    return(coverage)
}
