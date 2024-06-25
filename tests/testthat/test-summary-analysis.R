library(testthat)
source("../../R/utils.R")

# Test for summary_analysis function
test_that("summary_analysis function computes summary statistics correctly", {
  # Create a sample data frame for testing
  subset_x <- data.frame(
    Sex = c("M", "F", "M", "F", "M", NA, "F", "M"),
    Age = c("10-19", "20-29", "30-39", "40-49", NA, "10-19", "20-29", "30-39"),
    TM_1 = c("cLOF", "Ms", NA, "intronic", "cLOF", "cLOF", "Ms", "Ms"),
    TM_2 = c(NA, "cLOF", "Ms", NA, "cLOF", "Ms", "intronic", NA)
  )
  
  # Call the summary_analysis function
  result <- summary_analysis(subset_x)
  
  # Check if the result is a data frame
  expect_s3_class(result, "data.frame")
  
  # Check if the number of patients is computed correctly
  expect_equal(result[1, "Number of patients"], 8)
  
  # Check if sex information availability and distribution are computed correctly
  expect_equal(result[1, "Sex information available"], 7)
  expect_equal(result[1, "Number of male patients"], 4)
  expect_equal(result[1, "Number of female patients"], 3)
  expect_equal(round(result[1, "% of male"], digits = 2), 57.14)
  expect_equal(round(result[1, "% of female"], digits = 2), 42.86)
  
  # Check if age information availability and distribution are computed correctly
  expect_equal(result[1, "Age information available"], 7)
  expect_equal(result[1, "Number of age 10-19 years"], 2)
  expect_equal(result[1, "Number of age 20-29 years"], 2)
  expect_equal(result[1, "Number of age 30-39 years"], 2)
  # Add more checks for other age groups
  
  # Check if mutation type availability and distribution are computed correctly
  expect_equal(result[1, "Mutation type avaiable"], 12)
  expect_equal(result[1, "Number of cLOF"], 5)
  expect_equal(result[1, "Number of missense"], 5)
  expect_equal(result[1, "Number of intronic"], 2)
  expect_equal(round(result[1, "% of cLOF mutations"], digits = 2), 41.67)
  expect_equal(round(result[1, "% of missense mutations"], digits = 2), 41.67)
  expect_equal(round(result[1, "% of intronic mutations"], digits = 2), 16.67)
})
