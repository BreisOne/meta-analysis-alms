library(testthat)
source("../../R/utils.R")

# Define a test case for prevalence_statistics_chisq_global function
test_that("prevalence_statistics_chisq_global performs chi-square test correctly", {
  # Create sample data for testing
  update8forR <- data.frame(
                            VI = c(1, 0, 1, 1, 1, 0),
                            MT = c(0, 1, 0, 0, 0, 1),
                            HL = c(1, 1, 1, 0, 0, 0),
                            HRT = c(1, 1, 0, 0, 1, 1),
                            LIV = c(0, 1, 1, 0, 1, 1),
                            REN = c(1, 0, 1, 0, 1, 0),
                            PUL = c(1, 1, 1, 1, 0, 0),
                            SHS = c(1, 1, 1, 0, 1, 1),
                            REP = c(1, 0, 0, 1, 1, 1),
                            TYD = c(1, 1, 1, 0, 0, 1),
                            MEND = c(1, 1, 0, 1, 1, 1),
                            ABFING = c(1, 0, 1, 1, 1, 0),
                            INT = c(1, 1, 1, 1, 0, 1),
                            SCO = c(1, 1, 1, 0, 1, 0),
                            NER = c(1, 0, 1, 1, 1, 1),
                            ALO = c(1, 1, 0, 1, 0, 0)
                          )
    # Call the function
  result <- prevalence_statistics_chisq_global(update8forR)
  
  # Check if the result is a list
  expect_type(result, "double")
  
})