library(testthat)
source("../../R/utils.R")

test_that("prevalence_statistics_chisq_top returns only a numeric p-value", {
  # Create sample data for testing
  update8forR <- data.frame(
    VI = c(1, 0, 1, 1, 1, 0),
    MT = c(0, 1, 0, 0, 0, 1),
    HL = c(1, 1, 1, 0, 0, 0),
    HRT = c(1, 1, 0, 0, 1, 1),
    LIV = c(0, 1, 1, 0, 1, 1)
  )
  
  # Call the function
  result <- prevalence_statistics_chisq_top(update8forR)
  
  # Check if the result is a numeric value
  expect_type(result, "double")
})
