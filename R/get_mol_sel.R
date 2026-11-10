#' Get Model Selection Criteria
#'
#' Calculates model selection criteria such as AIC, BIC, DIC, and WAIC from a
#' Stan model fit.
#'
#' @param y A numeric vector representing the observed data.
#' @param mod_type Character string indicating the model type ("zi" or "h").
#' @param distri Character string specifying the distribution ("poi", "nb", "gp", or "cmp").
#' @param stan_fit A \code{stanfit} object returned by \code{get_stanfit}.
#' @param M CMP cutoff used in fitting. Retrieved from the fit returned by
#'   get_stanfit, or supplied explicitly for a fit made directly with Stan.
#' @return A summary table of model selection criteria, including:
#' \describe{
#'   \item{EAIC}{Expected Akaike Information Criterion (AIC).}
#'   \item{EBIC}{Expected Bayesian Information Criterion (BIC).}
#'   \item{DIC}{Deviance Information Criterion (DIC).}
#'   \item{WAIC1}{First version of Watanabe-Akaike Information
#'   Criterion (WAIC).}
#'   \item{WAIC2}{Second version of Watanabe-Akaike Information
#'   Criterion (WAIC).}
#' }
#' The criteria are returned in a one-row data frame.
#'
#' @examples
#' \donttest{
#'   # Generate toy data
#'   y_data <- data_simu(n = 60, alpha = 0.5, rho = 0.3, lambda = 5,
#'                       mod_type = "zi", distri = "poi")
#'
#'   # Fit a small Stan model (may take > 5s on first compile)
#'   stan_fit <- get_stanfit(mod_type = "zi", distri = "poi", y = y_data)
#'
#'   # Get model selection criteria
#'   get_mod_sel(y = y_data, mod_type = "zi", distri = "poi",
#'               stan_fit = stan_fit)
#' }
#'
#' @import rstan
#' @import knitr
#' @import matrixStats
#' @export
get_mod_sel <- function(y, mod_type, distri, stan_fit, M = NULL) {

  if (!is.numeric(y) || any(y != as.integer(y))) {
    stop("The parameter 'y' must be a numeric vector of integers.
         Please provide valid integer values.")
  }
  if (!(mod_type %in% c("zi", "h"))) {
    stop("The parameter 'mod_type' must be either 'zi' for zero-inflated
         or 'h' for hurdle.")
  }
  if (!(distri %in% c("poi", "nb", "gp", "cmp"))) {
    stop("The parameter 'distri' must be one of: 'poi', 'nb', 'gp', 'cmp'.")
  }
  if (mod_type == "h" && distri %in% c("gp", "cmp")) {
    stop(paste0("Hurdle models ('h') are not implemented for distribution '", distri,
                "'. Only Zero-Inflated ('zi') models are supported for GP and CMP."))
  }
  if (!inherits(stan_fit, "stanfit")) {
    stop("The parameter 'stan_fit' must be a valid 'stanfit' object
         returned by Stan.")
  }

  if (distri == "cmp") {
    all_criteria <- mod_sel_criteria(y, mod_type, distri, stan_fit, M = M)
    names(all_criteria) <- c("EAIC", "EBIC", "DIC", "WAIC1", "WAIC2")
    return(all_criteria)
  }

  aic = data.frame(extract(stan_fit, pars = "aic"))
  eaic = mean(aic[, 1])

  bic = data.frame(extract(stan_fit, pars = 'bic'))
  ebic = mean(bic[, 1])

  if (distri == 'poi') {
    para.hat = summary(stan_fit, pars = c("alpha", "rho", "lambda"))$summary
    alpha.hat = para.hat[1, 1]
    rho.hat = para.hat[2, 1]
    lambda.hat = para.hat[3, 1]

    logphat = sum(get_loglik(y, alpha = alpha.hat, rho = rho.hat,
                             lambda = lambda.hat, disp = NA,
                             mod_type, distri))
  } else {
    para.hat = summary(stan_fit,
                       pars = c("alpha", "rho", "lambda", "phi"))$summary
    alpha.hat = para.hat[1, 1]
    rho.hat = para.hat[2, 1]
    lambda.hat = para.hat[3, 1]
    disp.hat = para.hat[4, 1]

    logphat = sum(get_loglik(y, alpha = alpha.hat, rho = rho.hat,
                             lambda = lambda.hat, disp = disp.hat,
                             mod_type, distri))
  }

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

  return(all_criteria)
}
