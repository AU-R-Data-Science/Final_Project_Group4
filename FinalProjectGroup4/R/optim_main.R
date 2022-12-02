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

#' Estimate confidence intervals for the coefficients of logistic regression
#' @description This function bootstraps the input data to create sample populations to estimate the coefficients of logistic regression.
#' @param y A \code{double} value of the vector containing the response of interest.
#' @param X An \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors.
#' @param alpha A \code{double} value indicating the significance level of the confidence intervals.
#' @param rounds A \code{double} value indicating how many times to perform the bootstrap.
#' @return A \code{list} containing the following objects:
#' \describe{
#'  \item{boot_ci_int}{The estimated confidence interval for the coefficient of the intercept}
#'  \item{boot_ci_slope}{The estimated confidence interval for the coefficient of the slope}
#' }
#' @author Micheal Stewart Jackson
#' @export
bootstrap_ci <- function(y, X, alpha = .05, rounds = 20)
{
  beta_mat <- matrix(NA, rounds, 2)
  data <- data.frame("y" = y, "X" = X)
  for (b in 1:rounds)
  {
    boot_data <- data[sample(1:nrow(data), nrow(data), replace = TRUE), ]
    y2 <- boot_data[ , 1]
    X2 <- boot_data[ , 2]
    beta_mat[b, ] <- ls_optim(y2, X2)$beta_hat$par
  }
  boot_ci_int <- quantile(beta_mat[ , 1], c(alpha/2, 1-alpha/2))
  boot_ci_slope <- quantile(beta_mat[ , 2], c(alpha/2, 1-alpha/2))
  return(list("Bootstrap Confidence Interval For Intercept Coefficient"=boot_ci_int, "Bootstrap Confidence Interval For Slope Coefficient"=boot_ci_slope))
}

#' Create a confusion matrix and output its metrics
#' @description This function provides and uses a confusion matrix to assess how closely an optimization technique models the data presented.
#' @param y A \code{double} value of the vector containing the response of interest.
#' @param X An \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors.
#' @param alpha A \code{double} value indicating the significance level of the confidence intervals for the logistic regression coefficients.
#' @param cutoff A \code{double} value that sets the probability which determines the predicticted classification of the data set.
#' @return A \code{list} containing the following objects:
#' \describe{
#'  \item{Confusion Matrix}{The confusion matrix for the logistic regression}
#'  \item{Prevalence}{The percentage of positive cases in the data set}
#'  \item{Accuracy}{The percentage of cases predicted correctly by the model}
#'  \item{Sensitivity}{The percentage of positive cases predicted correctly by the model}
#'  \item{Specificity}{The percentage of negative cases predicted correctly by the model}
#'  \item{False Discovery Rate}{The percentage of cases predicted positive by the model that are actually negative}
#'  \item{Diagnostic Odds Ratio}{The ratio of the probability of a positive predicted as a negative to the probability of a negative predicted as a negative}
#' }
#' @author Micheal Stewart Jackson
#' @export
conf_matrix <- function(y, X, alpha = .05, cutoff = 0.5)
{


  int_coeff <- mean(bootstrap_ci(y, X, alpha)$"Bootstrap Confidence Interval For Intercept Coefficient")
  slope_coeff <- mean(bootstrap_ci(y, X, alpha)$"Bootstrap Confidence Interval For Slope Coefficient")

  log_odds <- int_coeff + slope_coeff * X
  probabilities <- exp(log_odds)/(1+exp(log_odds))
  predict <- ifelse(probabilities > cutoff, yes = 1, no = 0)

  tp <- 0 # Number of true positives
  fn <- 0 # Number of false negatives
  fp <- 0 # Number of false positives
  tn <- 0 # Number of true negatives

  for (i in 1:length(y))
  {
    if (y[i] < .5 & predict[i] < .5)
    {
      tn <- tn +1
    }
    if (y[i] < .5 & predict[i] > .5)
    {
      fp <- fp +1
    }
    if (y[i] > .5 & predict[i] > .5)
    {
      tp <- tp +1
    }
    if (y[i] > .5 & predict[i] < .5)
    {
      fn <- fn +1
    }
  }

  row1 <- c(tp, fn)
  row2 <- c(fp, tn)
  cm <-  rbind(row1, row2)# confusion matrix
  lab <- c(1,0)
  rownames(cm) <- lab
  colnames(cm) <- lab

  prev <- (tp+fp)/(tp+fp+fn+tn)
  acc <- (tp+tn)/(tp+fp+tn+fn)
  tpr <- tp/(tp+fn)
  tnr <- tn/(tn+fp)
  fdr <- fp/(fp+tp)
  dor <- (tpr/(1-tnr))/((1-tpr)/tnr)


  return(list("Confusion Matrix" = cm, "Prevalence" = prev, "Accuracy" = acc, "Sensitivity" = tpr, "Specificity" = tnr, "False Discovery Rate" = fdr, "Diagnostic Odds Ratio" = dor))
}


#' Create a plot of accuracy from a confusion matrix
#' @description This function shows plot of accuracy from a confusion matrix evaluated over a grid of cut-off values for prediction going from 0.1 to 0.9
#' @param y A \code{double} value of the vector containing the response of interest.
#' @param X An \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors.
#' @param alpha A \code{double} value indicating the significance level of the confidence intervals for the logistic regression coefficients.
#' @return A \code{list} containing the following objects:
#' \describe{
#'  \item{model}{}
#' }
#' @author Kayla Gallman
#' @export
plot_metrics <- function(y, X, alpha = 0.5)
{

  cut_off <- seq(0.1, 0.9, by = 0.1)
  for (i in 1:length(cut_off)){
    co <- cut_off[i]
    Confusion <- conf_matrix(y,X, alpha = 0.5, cutoff = co)
  }
  plot(cut_off, Confusion$Accuracy)

}

set.seed(1)
y <- sample(c(0,1), size = 100, replace = TRUE)
X <- round(runif(100, 18, 80))


#' Create a fitted logistic curve
#' @description This function shows a visual of the logistic regression
#' @param y A \code{double} value of the vector containing the response of interest.
#' @param X An \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors.
#' @param beta A \code{double} value that has previously been found with the ls_optim function.
#' @return A \code{list} containing the following objects:
#' \describe{
#'  \item{model}{fit of the logistic regression model}
#' }
#' @author Kayla Gallman
#' @export
log_curve<- function(X, y, beta)
{

  #fit logistic regression model
  model <- ls_optim(y,X)

  #define new data frame that contains predictor variable
  newdata <- data.frame(X=seq(min(X), max(X),len=500))

  #use fitted model to predict values of vs

  p <- 1/exp(-y*X)

  #plot logistic regression curve
  plot(y ~ X, col="steelblue")
  lines(y ~ X, newdata, lwd=2)

}
