# Define a test case for phenotype_analysis function
test_that("phenotype_analysis calculates statistics correctly", {
  # Create sample data frame for testing
  df <- data.frame(VI = c(1, 0, NA, 1, 1, 0),
                   LIV = c(0, 1, NA, 1, 0, 0))
  
  # Call the function
  result_VI <- phenotype_analysis(df$VI)
  result_LIV <- phenotype_analysis(df$LIV)
  
  
  # Check if the result is a data frame
  expect_s3_class(result_VI, "data.frame")
  expect_s3_class(result_LIV, "data.frame")
  
  # Check if the result has expected column names
  expect_equal(names(result_VI), c("count", "yes", "no", "na", "percent", "upper", "lower"))
  
  # Check if the count of non-missing values is correct
  expect_equal(result_VI$count, 5)
  
  # Check if the count of "yes" values for visual impairment is correct
  expect_equal(result_LIV$yes, 2)
  expect_equal(result_VI$yes, 3)
  
  # Check if the count of "no" values for visual impairment is correct
  expect_equal(result_VI$no, 2)
  
  # Check if the count of missing values for visual impairment is correct
  expect_equal(result_LIV$na, 1)
  
  # Check if the percentage of "yes" values for visual impairment is correct
  expect_equal(result_LIV$percent, 40)
  
  # Check if the upper confidence interval for visual impairment is calculated correctly
  expect_equal(round(result_VI$upper, digits = 2), 88.19)
  
  # Check if the lower confidence interval for visual impairment is calculated correctly
  expect_equal(result_LIV$lower, 11.81)

})