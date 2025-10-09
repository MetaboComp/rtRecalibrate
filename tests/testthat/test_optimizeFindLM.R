test_that("optimizeFindLM", {

  #Check that not both PT_filepath and PT submitted by user
  expect_error(optimizeFindLM(PT_filepath="Ö:/test.rds",
                               PT = peak_table,
                               split='@',
                               minrttocheck=40,
                               Prefilterintensity=7000,
                               minLM=100,
                               maxLM=200))
  
  #Check that PT_filepath file exists
  expect_error(optimizeFindLM(PT_filepath="Ö:/test.rds",
                               split='@',
                               minrttocheck=40,
                               Prefilterintensity=7000,
                               minLM=100,
                               maxLM=200))
  
  #Check number of features is higher than minLM
  expect_error(optimizeFindLM(PT=peak_table[,seq_len(50)],
                              split='@',
                              minrttocheck=40,
                              Prefilterintensity=7000,
                              minLM=100,
                              maxLM=200))
  
  #Check that missing Prefilterintensity throws error
  expect_error(optimizeFindLM(PT=peak_table,
                              split='@',
                              minrttocheck=40,
                              Prefilterintensity,
                              minLM=100,
                              maxLM=200))
  
})