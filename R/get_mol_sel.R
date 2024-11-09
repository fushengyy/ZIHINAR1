#' Get Model Selection Criteria
#'
#' Calculates model selection criteria such as AIC, BIC, DIC, and WAIC from a Stan model fit.
#'
#' @param y A numeric vector representing the observed data.
#' @param mod_type Character string indicating the model type ("zi" or "h").
#' @param distri Character string specifying the distribution ("poi" or "nb").
#' @param stan_fit A \code{stanfit} object returned by \code{get_stanfit}.
#' @return A summary table of model selection criteria, including:
#' \describe{
#'   \item{EAIC}{Expected Akaike Information Criterion (AIC).}
#'   \item{EBIC}{Expected Bayesian Information Criterion (BIC).}
#'   \item{DIC}{Deviance Information Criterion (DIC).}
#'   \item{WAIC1}{First version of Watanabe-Akaike Information Criterion (WAIC).}
#'   \item{WAIC2}{Second version of Watanabe-Akaike Information Criterion (WAIC).}
#' }
#' The summary is printed in a table format for easy interpretation.
#'
#' @examples
#' \dontrun{
#'   # Generate simulated data
#'   y_data <- data_simu(n = 100, alpha = 0.5, rho = 0.3, theta = c(5), mod_type = "zi", distri = "poi")
#'
#'   # Fit the model using Stan
#'   stan_fit <- get_stanfit(mod_type = "zi", distri = "poi", y = y_data)
#'
#'   # Get model selection criteria
#'   get_mod_sel(y = y_data, mod_type = "zi", distri = "poi", stan_fit = stan_fit)
#' }
#'
#' @import knitr
#' @export
get_mod_sel <- function(y, mod_type, distri, stan_fit) {

  install_if_missing("rstan")
  install_if_missing("knitr")

  if (!is.numeric(y) || any(y != as.integer(y))) {
    stop("The parameter 'y' must be a numeric vector of integers. Please provide valid integer values.")
  }
  if (!(mod_type %in% c("zi", "h"))) {
    stop("The parameter 'mod_type' must be either 'zi' for zero-inflated or 'h' for hurdle.")
  }
  if (!(distri %in% c("poi", "nb"))) {
    stop("The parameter 'distri' must be either 'poi' for Poisson or 'nb' for Negative Binomial.")
  }
  if (!inherits(stan_fit, "stanfit")) {
    stop("The parameter 'stan_fit' must be a valid 'stanfit' object returned by Stan.")
  }

  aic = data.frame(extract(stan_fit, pars = "aic"))
  eaic = mean(aic[, 1])

  bic = data.frame(extract(stan_fit, pars = 'bic'))
  ebic = mean(bic[, 1])

  if (distri == 'poi') {
    para.hat = summary(stan_fit, pars = c("alpha", "rho", "lambda"))$summary
    alpha.hat = para.hat[1, 1]
    rho.hat = para.hat[2, 1]
    theta.hat = para.hat[3, 1]
  } else {
    para.hat = summary(stan_fit, pars = c("alpha", "rho", "lambda", "phi"))$summary
    alpha.hat = para.hat[1, 1]
    rho.hat = para.hat[2, 1]
    theta.hat = c(para.hat[3, 1], para.hat[4, 1])
  }

  logphat = sum(get_loglik(y, alpha = alpha.hat, rho = rho.hat, theta = theta.hat, mod_type, distri))
  ll = data.frame(extract(stan_fit, pars = "ll"))
  pdic = 2 * (logphat - mean(ll[, 1]))
  dic = -2 * logphat + 2 * pdic

  lik = data.frame(extract(stan_fit, pars = "lik"))
  loglik = data.frame(extract(stan_fit, pars = "log_lik"))
  pwaic1 = 2 * sum(log(colMeans(lik)) - colMeans(loglik))
  pwaic2 = sum(colVars(as.matrix(loglik)))
  lppd = sum(log(colMeans(lik)))
  waic1 = -2 * (lppd - pwaic1)
  waic2 = -2 * (lppd - pwaic2)

  all_criteria = data.frame(eaic, ebic, dic, waic1, waic2)
  colnames(all_criteria) = c("EAIC", "EBIC", "DIC", "WAIC1", "WAIC2")
  print(knitr::kable(all_criteria, format = "simple", digits = 4))
}
