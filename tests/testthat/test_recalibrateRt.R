test_that("rtRecalibrate",{
  
  #Checking that not choosing a method throws an error
  expect_error(suppressWarnings(
    MsExpObj_result <- rtRecalibrate:::rtRecalibrate(MsExpObj,
                             landmarks_final,
                             method = c('GAM', 'loess'),
                             roiExpandRt = 15,
                             roiPpm = 40,
                             toleranceRt = 15,
                             ppm = 5,
                             tolerance = 0,
                             peakwidth = c(7, 50),
                             bs = "tp",
                             span = 0.5,
                             ssqRatio = 3,
                             zeroWeight = 10,
                             BPPARAM = bpparam())))
  
  #Checking that choosing wrong method throws error
  expect_error(rtRecalibrate:::rtRecalibrate(MsExpObj,
                             landmarks_final,
                             method = c('wrong'),
                             roiExpandRt = 15,
                             roiPpm = 40,
                             toleranceRt = 15,
                             ppm = 5,
                             tolerance = 0,
                             peakwidth = c(7, 50),
                             bs = "tp",
                             span = 0.5,
                             ssqRatio = 3,
                             zeroWeight = 10,
                             BPPARAM = bpparam()))
  
  #Checking that processing an MSnExp object works
  expect_no_error(suppressWarnings(
    MSnExpObj_result <- rtRecalibrate:::rtRecalibrate(MSnExpObj,
                               landmarks_final,
                               method = c('GAM'),
                               roiExpandRt = 15,
                               roiPpm = 40,
                               toleranceRt = 15,
                               ppm = 5,
                               tolerance = 0,
                               peakwidth = c(7, 50),
                               bs = "tp",
                               span = 0.5,
                               ssqRatio = 3,
                               zeroWeight = 10,
                               plot=FALSE,
                               BPPARAM = bpparam())))
  
  #Checking that processing an MsExperiment object works
  expect_no_error(suppressWarnings(
    MsExpObj_result <- rtRecalibrate:::rtRecalibrate(MSnExpObj,
                               landmarks_final,
                               method = c('GAM'),
                               roiExpandRt = 15,
                               roiPpm = 40,
                               toleranceRt = 15,
                               ppm = 5,
                               tolerance = 0,
                               peakwidth = c(7, 50),
                               bs = "tp",
                               span = 0.5,
                               ssqRatio = 3,
                               zeroWeight = 10,
                               BPPARAM = bpparam())))
  
  
  
})