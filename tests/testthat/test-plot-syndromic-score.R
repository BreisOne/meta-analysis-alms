library(ggplot2)
library(testthat)
source("../../R/utils.R")

test_that("plot_syndromic_score generates a plot object", {
  # Create sample data for testing
  df <- data.frame(SS = rnorm(100))
  
  # Call the function
  plot <- plot_syndromic_score(df, "Frequency")
  
  # Check if the result is a ggplot object
  expect_s3_class(plot, "ggplot")
})
