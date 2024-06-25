library(testthat)
library(dplyr)
source("../../R/utils.R")

# Create sample data for testing
update8forR <- update8forR <- data.frame(
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

# Define a test case for prevalence_symptons_global function
test_that("prevalence_symptons_global calculates prevalence correctly", {
  # Call the function
  result <- prevalence_symptons_global(update8forR)
  
  # Check if the result is a data frame
  expect_s3_class(result, "data.frame")
  
  # Check if the result has expected number of rows (phenotypes)
  expect_equal(nrow(result), 16)
  
  # Check if the result has expected number of columns (1 for percentage)
  expect_equal(ncol(result), 1)
  
  # Check if the result is sorted in descending order of prevalence percentages
  expect_true(all(diff(result$Percent) <= 0))
})
