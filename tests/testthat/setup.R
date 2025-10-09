set.seed(2024)

################# Find landmarks #################
landmark_suggestions <- optimizeFindLM(PT=peak_table,
                                       split='@',
                                       Prefilterintensity=7000,
                                       minLM = 40)

landmarks_final <- landmark_suggestions[[1]][[2]]

################# Recalibrate rt of .mzML files #################
MsBackendMetaboLights::mtbls_sync_data_files(mtblsId = "MTBLS8735",
                      assayName = "a_MTBLS8735_LC-MS_positive_hilic_metabolite_profiling.txt")

mzMLFiles <- MsBackendMetaboLights::mtbls_cached_data_files()$rpath
mzMLFiles_QC <- mzMLFiles[grepl("QC", mzMLFiles)]

MsExpObj <- MsExperiment::readMsExperiment(mzMLFiles_QC)
MSnExpObj <- MSnbase::readMSData(mzMLFiles_QC,
                                 mode = "onDisk")