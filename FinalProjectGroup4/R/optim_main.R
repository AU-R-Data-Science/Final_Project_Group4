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

#' Estimate linear model through numerical optimization
#' @description This function provides a least-squares regression through the minimization of beta values in a specified loss function.
#' @param y A \code{double} value of the vector containing the response of interest.
#' @param X An \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors.
#' @return A \code{list} containing the following objects:
#' \describe{
#'  \item{beta_hat}{The estimated coefficients of the linear regression}
#'  \item{y}{The \code{double} vector containing the response used for the estimation}
#'  \item{X}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation}
#' }
#' @author Mauren Baker
#' @export
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
