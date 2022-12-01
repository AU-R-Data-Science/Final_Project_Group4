library(boot)

ls_obj <- function(y, X)
{
  y <- as.vector(y)
  X <- as.matrix(X)
  solve(t(X)%*%X)%*%t(X)%*%y
}

loss <- function(y, X, beta)
{
  p <- c()
  bh <- c()
  n <- length(y)

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
ls_optim <- function(y, X)
{
  n <- length(y)
  X <- cbind(X, c(rep(1, n)))
  beta_est <- optim(par = ls_obj(y, X), loss, y = y, X = X)

  out <- list("beta_hat" = beta_est, "response" = y, "predictors" = X)
  return(out)
}


bootstrap_ci <- function(alpha, rounds = 20)
{
  beta_mat <- matrix(NA, B, 2)
    for (b in 1:rounds)
    {
     boot_data <- data[sample(1:nrows(data), nrows(data), replace = TRUE), ]
     beta_mat[b, ] <- lm(y~x, data = boot_data)$coefficients
    }

  return(list(bootstrap_ci = boot_ci))
}

#' Create a confusion matrix and output its metrics
#' @description This function provides and uses a confusion matrix to assess how closely an optimization technique models the data presented.
#' @param y A \code{double} value of the vector containing the response of interest.
#' @param X An \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors.
#' @param cutoff A \code{double} value that sets the probability which determines the predicticted classification of the data set.
#' @return A \code{list} containing the following objects:
#' \describe{
#'  \item{Prevalence}{The percentage of positive cases in the data set}
#'  \item{Accuracy}{The percentage of cases predicted correctly by the model}
#'  \item{Sensitivity}{The percentage of positive cases predicted correctly by the model}
#'  \item{Specificity}{The percentage of negative cases predicted correctly by the model}
#'  \item{False Discovery Rate}{The percentage of cases predicted positive by the model that are actually negative}
#'  \item{Diagnostic Odds Ratio}{The ratio of the probability of a positive predicted as a negative to the probability of a negative predicted as a negative}
#' }
#' @author Micheal Stewart Jackson
#' @export
  conf_matrix <- function(y, X, cutoff = 0.5)
  {
  ls_optim(y, X)
    predict <- ifelse(beta_est > cutoff, yes = 1, no = 0)
  cm <- table(y, predict) # confusion matrix

  tp <- cm[1, 1] # Number of true positives
  fn <- cm[1, 2] # Number of false negatives
  fp <- cm[2, 1] # Number of false positives
  tn <- cm[2, 2] # Number of true negatives

  prev <- (tp+fp)/(tp+fp+fn+tn)
  acc <- (tp+tn)/(tp+fp+tn+fn)
  tpr <- tp/(tp+fn)
  tnr <- tn/(tn+fp)
  fdr <- fp/(fp+tp)
  dor <- (tpr/(1-tnr))/((1-tpr)/tnr)


return(list("Prevalence" = prev, "Accuracy" = acc, "Sensitivity" = tpr, "Specificity" = tnr, "False Discovery Rate" = fdr, "Diagnostic Odds Ratio" = dor))
  }


  log_curve <- function(y, X)
  {

    #fit logistic regression model
    model <- ls_optim(y,X)

    #define new data frame that contains predictor variable
    lr <- data.frame(X=seq(min(X), max(X),len=100))

    #plot logistic regression curve
    plot(y ~ seq(X), col="steelblue")
    lines(y ~ X, lr, lwd=2)

  }

  plot_metrics <- function()
  {

  }
