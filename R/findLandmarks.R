#' findLandmarks - A function for finding landmarks within
#'  a resulting data frame from xcms
#'
#' @param dat Data frame with features in columns and samples in rows
#' @param mzrtdf Data frame with mz values of all features in column 1
#'  and retention time values of all features in column 2 named mzvec and rtvec
#' @param mzdiff How large whitespace in mz a feature needs, typically 0.1
#' @param rtdiff How large whitespace in rt a feature needs,
#'  can vary quite a lot, but ~20 is sufficient
#' @param minIntensity Value indicating what min intensity is needed
#'  for a feature to qualify as an LM
#' @param minrttocheck Retention time value where noisy features
#'  (that usually comes in the beginning of a run) can be ignored
#' @param allowedmissingness Double between 0 - 1 specifying how 
#' much missingness allowed for LaMas within all samples
#'
#' @return Returning potential Landmarks
#' 

findLandmarks <- function(dat,
                          mzrtdf,
                          mzdiff,
                          rtdiff,
                          minIntensity,
                          minrttocheck,
                          allowedmissingness)
{
  Landmarks <- c(NA,NA)
  #Calculate missingness in each feature
  PropNAsinData <- apply(dat, 2, function (x) sum(is.na(x))/length(x))

  #Only check features with <20% missingness
  FeaturesWithFewNAs <- dat[,PropNAsinData < allowedmissingness]
  mzrtwithfewNAs <- mzrtdf[PropNAsinData < allowedmissingness, ]
  Intensityvector <- apply(FeaturesWithFewNAs, 2, mean, na.rm=TRUE)

  LMcandidates <- FeaturesWithFewNAs[, Intensityvector > minIntensity]
  mzrtcandidates <- mzrtwithfewNAs[Intensityvector > minIntensity,]
  LMcandidates <- LMcandidates[, mzrtcandidates$rtvec > minrttocheck]
  mzrtcandidates <- mzrtcandidates[mzrtcandidates$rtvec > minrttocheck, ]

  if(ncol(LMcandidates)>10){
    Landmarkmzs <- c()
    Landmarkrts <- c()

    for(i in seq_len(ncol(LMcandidates))){
      mzofInterest <- mzrtcandidates$mzvec[i]
      rtofInterest <- mzrtcandidates$rtvec[i]
      Featsinregion <- sum(mzrtdf$mzvec>mzofInterest-mzdiff &
                           mzrtdf$mzvec<mzofInterest+mzdiff &
                           mzrtdf$rtvec>rtofInterest-rtdiff &
                           mzrtdf$rtvec<rtofInterest+rtdiff)

      if(Featsinregion==1){
        Landmarkmzs <- c(Landmarkmzs, mzofInterest)
        Landmarkrts <- c(Landmarkrts, rtofInterest)
      }
    }

    Landmarks <- data.frame(Landmarkmzs, Landmarkrts)
    colnames(Landmarks) <- c("MZ", "RT")
  }

return(Landmarks)
}
