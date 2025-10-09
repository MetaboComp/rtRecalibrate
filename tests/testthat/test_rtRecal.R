test_that("rtRecal", {
  # mtbls_sync_data_files(mtblsId = "MTBLS8735",
  #     assayName = "a_MTBLS8735_LC-MS_positive_hilic_metabolite_profiling.txt")
  # 
  # mzMLFiles <- mtbls_cached_data_files()$rpath
  
  #Checking that non-existing file/directory throws and error
  expect_error(rtRecal("Ö:/Asdf.fake",
                       landmarks_final,
                       method = c('GAM'),
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
                       save = TRUE))
  
  #Checking that non-mzML file directory throws an error
  expect_error(rtRecal(paste0(.libPaths(),"/rtRecalibrate/landmarks_final.rda"),
                       landmarks_final,
                       method = c('GAM'),
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
                       save = TRUE))
  
  #Checking that not choosing a method throws an error
  expect_error(rtRecal(mzMLFiles[1],
                       landmarks_final,
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
                       save = TRUE))
  
  #Checking that wrong format landmarks throws an error
  landmarks_bad <- as.data.frame(matrix(ncol=3, nrow=10))
  colnames(landmarks_bad) <- c("MZ", "RT")
  
  expect_error(rtRecal(mzMLFiles[1],
                       landmarks_bad,
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
                       save = TRUE))
  
  landmarks_bad <- matrix(ncol=2, nrow=10)
  colnames(landmarks_bad) <- c("MZ", "RT")
  
  expect_error(rtRecal(mzMLFiles[1],
                       landmarks_bad,
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
                       save = TRUE))
})