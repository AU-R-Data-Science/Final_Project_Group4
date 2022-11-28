least_squares <- function(response, predictors)
{
  response <- as.vector(response)
  predictors <- as.matrix(predictors)

  inv_pred <- solve(t(predictors) %*% predictors)
  beta <- inv_pred %*% t(predictors) %*% response
  return(beta)
}

optimize <- function(response, predictors, beta)
{
  for (i in predictors)
    {
      p[i] <- 1 / (1 + exp(t(-1 * predictors[i]) * beta))
    }
}

bootstrap_ci <- function()
{
  
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
