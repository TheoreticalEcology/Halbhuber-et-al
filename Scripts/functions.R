calculate_r_squared <- function(true_values, predictions) {
  res_ss <- sum((true_values - predictions)^2)
  total_ss <- sum((true_values - mean(true_values))^2)
  R_squared <- 1 - res_ss / total_ss
  return(R_squared)
}