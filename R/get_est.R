#' Get Parameter Estimates from Stan Model Fit
#'
#' Extracts parameter estimates from a Stan model fit, including mean, median,
#' standard deviation, and HPD intervals.
#'
#' @param distri Character string specifying the distribution. Options
#' are "poi" or "nb".
#' @param stan_fit A \code{stanfit} object returned by \code{get_stanfit}.
#' @return A summary of the parameter estimates.
#'
#' @examples
#' \donttest{
#'   # Generate toy data
#'   y_data <- data_simu(n = 60, alpha = 0.5, rho = 0.3, theta = c(5),
#'                       mod_type = "zi", distri = "poi")
#'
#'   # Fit a small Stan model (may take > 5s on first compile)
#'   stan_fit <- get_stanfit(mod_type = "zi", distri = "poi", y = y_data)
#'
#'   # Get parameter estimates from the Stan model fit
#'   get_est(distri = "poi", stan_fit = stan_fit)
#' }
#'
#' @import knitr
#' @import rstan
#' @importFrom coda HPDinterval as.mcmc as.mcmc.list
#' @export
get_est <- function(distri, stan_fit) {

  if (!inherits(stan_fit, "stanfit")) {
    stop("The parameter 'stan_fit' must be a valid 'stanfit' object
         returned by Stan.")
  }
  if (!(distri %in% c("poi", "nb"))) {
    stop("The parameter 'distri' must be either 'poi' for Poisson
         or 'nb' for Negative Binomial.")
  }

  if (distri == "poi") {
    qoi <- c("alpha", "rho", "lambda")
  } else {
    qoi <- c("alpha", "rho", "lambda", "phi")
  }

  est_all <- summary(stan_fit, pars = qoi)$summary
  selected_columns <- c("mean", "sd", "2.5%", "50%", "97.5%", "Rhat")
  est_selected <- est_all[, selected_columns]

  colnames(est_selected) <- c("Mean", "SD", "Q2.5", "Median", "Q97.5", "Rhat")
  est_selected <- est_selected[, c("Mean", "SD", "Median",
                                   "Q2.5", "Q97.5", "Rhat")]

  param_samples <- rstan::extract(stan_fit, pars = qoi)

  if (is.list(param_samples) && length(param_samples) > 1) {
    param_samples_df <- do.call(cbind, lapply(param_samples, function(chain) {
      as.data.frame(chain)
    }))

    combined_mcmc <- coda::as.mcmc(param_samples_df)

    combined_mcmc_list <- coda::as.mcmc.list(list(combined_mcmc))

    hpd_intervals <- coda::HPDinterval(combined_mcmc_list, prob = 0.95)
  } else {
    stop("The extracted samples do not have the expected structure.")
  }

  hpd_df <- as.data.frame(hpd_intervals)
  colnames(hpd_df) <- c("95%_HPD_Lower", "95%_HPD_Upper")
  rownames(hpd_df) <- qoi

  est_result <- cbind(est_selected, hpd_df)

  return(est_result)
}
