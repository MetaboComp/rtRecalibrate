test_that("rtRecalibratemzML",{
  
  #Checking that non-existing file/directory throws and error
  expect_error(rtRecalibrate(c("Ö:/Asdf.fake",
                               mzMLFiles[1]),
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
  expect_error(rtRecal(c(paste0(.libPaths(),
                                "/rtRecalibrate/landmarks_final.rda"),
                         mzMLFiles[1]),
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

