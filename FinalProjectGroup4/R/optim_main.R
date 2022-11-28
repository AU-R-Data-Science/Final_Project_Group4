least_squares <- function(response, predictors, alpha)
{
  response <- as.vector(response)
  predictors <- as.matrix(predictors)

  n <- length(response)
  p <- dim(predictors)[2]
  df <- n - p

  inv_pred <- solve(t(predictors) %*% predictors)
  beta <- inv_pred %*% t(predictors) %*% response

  residual <- response - predictors %*% as.matrix(beta)
  sig2 <- (1/df) * t(residual) %*% residual
  var_b <- sig2 * inv_pred

  quantile <- 1 - alpha/2
  ci <- c(beta - qnorm(p = quantile) * sqrt(var_b), beta + qnorm(p = quantile) * sqrt(var_b))

  return(list(beta = beta, sigma2 = sig2, variance_beta = var_b, confidence_interval = ci))
}

optimize <- function(response, predictors, beta)
{
  for (i in predictors)
    {
      p[i] <- 1 / (1 + exp(t(-1 * predictors[i]) * beta))
    }
}

conf_matrix <- function()
{

}

log_curve <- function()
{

}

plot_metrics <- function()
{

}
