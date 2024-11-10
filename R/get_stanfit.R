#' Fit ZI-INAR(1) or H-INAR(1) Model using Stan
#'
#' This function fits a Zero-Inflated INAR(1) (ZI-INAR(1)) or Hurdle INAR(1)
#' (H-INAR(1)) model using Stan and returns the model fit.
#'
#' @param mod_type Character string indicating the model type.
#' Use "zi" for zero-inflated
#' models and "h" for hurdle models.
#' @param distri Character string specifying the distribution.
#' Options are "poi" for Poisson or "nb" for Negative Binomial.
#' @param y A numeric vector of integers representing the observed data.
#' @param n_pred Integer specifying the number of time points for future
#' predictions (default is 4).
#' @param thin Integer indicating the thinning interval for Stan sampling
#' (default is 2).
#' @param chains Integer specifying the number of Markov chains to run
#' (default is 1).
#' @param iter Integer specifying the total number of iterations per chain
#' (default is 2000).
#' @param warmup Integer specifying the number of warmup iterations per chain
#' (default is \code{iter/2}).
#' @param seed Numeric seed for reproducibility (default is NA).
#'
#' @return A \code{stanfit} object containing the Stan model fit.
#'
#' @examples
#' \dontrun{
#'   # Generate simulated data
#'   y_data <- data_simu(n = 100, alpha = 0.5, rho = 0.3, theta = c(5),
#'                       mod_type = "zi", distri = "poi")
#'
#'   # Fit the model using Stan
#'   stan_fit <- get_stanfit(mod_type = "zi", distri = "poi", y = y_data)
#'   print(stan_fit)
#' }
#'
#' @import rstan
#' @export
get_stanfit <- function(mod_type, distri, y, n_pred = 4,
                        thin = 2, chains = 1, iter = 2000, warmup = iter/2,
                        seed = NA) {

  install_if_missing("rstan")

  if (!is.numeric(y) || any(y != as.integer(y))) {
    stop("The parameter 'y' must be a numeric vector of integers.
         Please provide valid integer values.")
  }
  if (!(mod_type %in% c("zi", "h"))) {
    stop("The parameter 'mod_type' must be either 'zi' for zero-inflated
         or 'h' for hurdle.")
  }
  if (!(distri %in% c("poi", "nb"))) {
    stop("The parameter 'distri' must be either 'poi' for Poisson
         or 'nb' for Negative Binomial.")
  }


  stan_file <- if (dir.exists("stan")) "stan" else file.path("inst", "stan")
  if (mod_type == "zi" && distri == "poi") {
    stan_file <- file.path(stan_file, "ZIP-INAR1.stan")
  } else if (mod_type == "zi" && distri == "nb") {
    stan_file <- file.path(stan_file, "ZINB-INAR1.stan")
  } else if (mod_type == "h" && distri == "poi") {
    stan_file <- file.path(stan_file, "HP-INAR1.stan")
  } else {
    stan_file <- file.path(stan_file, "HNB-INAR1.stan")
  }

  stan_fit <- rstan::stan(file = stan_file,
                          data = list(y = c(y), T = length(y), ff = n_pred),
                          thin = thin,
                          chains = chains,
                          iter = iter,
                          warmup = warmup,
                          seed = seed)

  return(stan_fit)
}
