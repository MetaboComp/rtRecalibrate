#' Landmarks Final Dataset
#'
#' Landmarks found through use of the 'optimizeLM()' function performed on a
#' peak table derived from xcms preprocessing up until, but not including,
#' gap filling
#'
#' @name landmarks_final
#' @docType data
#' @format A data frame containing mz and rt for all landmarks where each row is
#' a landmark
#' \describe{
#'   \item{MZ}{Mass-to-charge ratio of landmark features}
#'   \item{RT}{Retention time of landmark features}
#'   ...
#' }
#' @usage data(landmarks_final)
"landmarks_final"

#' Peak Table Dataset
#'
#' Peak table generated in xcms using IPO and IPO2 to find optimized parameters
#' for xmcs functions: findChromPeaks, adjustRtime, groupChromPeaks and finally
#' featureValues and featureDefinitions to build the peak table
#'
#' @name peak_table
#' @docType data
#' @format A data frame peak table with sample in rows and features in columns
#' where mz and rt
#' @usage data(peak_table)
"peak_table"
