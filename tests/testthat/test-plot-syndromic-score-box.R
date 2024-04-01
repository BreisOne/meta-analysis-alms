library(ggplot2)

test_that("plot_syndromic_score_box generates a plot object", {
  # Create sample data for testing
  df <- data.frame(
    x_var = rep(c("Group A", "Group B"), each = 50),
    SS = c(rnorm(50), rnorm(50, mean = 0.5))
  )
  
  # Call the function
  plot <- plot_syndromic_score_box(df, "x_var")
  
  # Check if the result is a ggplot object
  expect_is(plot, "ggplot")
})
