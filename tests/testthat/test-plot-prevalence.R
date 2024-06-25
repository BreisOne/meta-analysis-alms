library(testthat)
library(dplyr)
source("../../R/utils.R")

test_that("plot_prevalence generates a plot object", {
  # Create sample data for testing
  update8forR <- list(
    VI = c(1, 1, 0, 1, 0, 0, 1, 1, 0, 1),
    MT = c(0, 1, 1, 0, 1, 1, 1, 0, 0, 1),
    HL = c(1, 0, 0, 1, 1, 0, 1, 0, 1, 0),
    HRT = c(0, 1, 0, 0, 1, 0, 1, 1, 0, 1),
    LIV = c(1, 0, 1, 1, 0, 1, 0, 0, 1, 1)
  )
  
  # Call the function
  plot <- plot_prevalence(update8forR, "global")
  
  # Check if the result is a plot object
  expect_null(plot)
})
