library(boot)

ls_obj <- function(y, X)
{
  solve(t(X)%*%X)%*%t(X)%*%y
}

loss <- function(y, X, beta)
{
  p <- c()
  bh <- c()
  n = length(y)
  #X <- cbind(X, c(rep(1, n)))
  for (i in 1:n)
  {
    p[i] <- 1 / (1 + exp(t(-X[i, ]) %*% beta))
  }

  for (i in 1:n)
  {
    bh[i] <- (-y[i] * log(p[i]) - (1 - y[i]) * log(1 - p[i]))
  }

  return(sum(bh))
}

optimize <- function(y, X)
{
  beta_est <- optim(par = ls_obj(y, X), loss, y = y, X = X)

  out <- list("beta_hat" = beta_est, "response" = y, "predictors" = X)
  return(out)
}

bootstrap_ci <- function(alpha, rounds = 20)
{
  beta_mat <- matrix(NA, B, )
  #for (b in 1:rounds)
  #{
   # boot_data <-
  # beta_mat[b, ] <- lm()$coefficients
  #}

 }

  return(list(bootstrap_ci = boot_ci))
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
