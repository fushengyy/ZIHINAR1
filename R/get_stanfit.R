#' Fit ZI-INAR(1) or H-INAR(1) Model using Stan
#'
#' This function fits a Zero-Inflated INAR(1) (ZI-INAR(1)) or Hurdle INAR(1)
#' (H-INAR(1)) model using Stan and returns the model fit.
#'
#' @param mod_type Character string indicating the model type.
#' Use "zi" for zero-inflated models and "h" for hurdle models.
#' @param distri Character string specifying the distribution.
#' Options are "poi" for Poisson, "nb" for Negative Binomial,
#' "gp" for Generalized Poisson, or "cmp" for Conway-Maxwell-Poisson.
#' GP and CMP innovations are supported only for zero-inflated models.
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
#' @param M CMP innovation cutoff. Defaults to the larger of 300 and
#' \code{max(y)}.
#' @param cmp_version CMP Stan implementation, either "fast" or "reference".
#'   Both use the same likelihood and priors. Defaults to "fast".
#' @param cmp_tail_tol Tolerance for the CMP omitted probability mass bound.
#'   Defaults to 1e-6.
#' @param tail_action Action when the CMP tail bound exceeds the tolerance.
#'   Use "warn" (default) or "error". Increase M and refit if flagged.
#' @param ... Additional arguments passed to \code{rstan::stan}.
#'
#' @return A \code{stanfit} object containing the Stan model fit.
#'
#' @examples
#' \donttest{
#'   # Generate toy data
#'   y_data <- data_simu(n = 60, alpha = 0.5, rho = 0.3, lambda = 5,
#'                       mod_type = "zi", distri = "poi")
#'
#'   # Fit the model using Stan (small config)
#'   stan_fit <- get_stanfit(mod_type = "zi", distri = "poi", y = y_data)
#'   print(stan_fit)
#' }
#'
#' \donttest{
#'   # Simulate and fit a zero-inflated CMP-INAR(1) series
#'   y_cmp <- data_simu(n = 60, alpha = 0.5, rho = 0.3, lambda = 2.5,
#'                      disp = 1.2, mod_type = "zi", distri = "cmp")
#'   fit_cmp <- get_stanfit(mod_type = "zi", distri = "cmp", y = y_cmp,
#'                          M = 300, cmp_version = "fast")
#'   print(fit_cmp)
#' }
#'
#' @import rstan
#' @export
get_stanfit <- function(mod_type, distri, y, n_pred = 4,
                        thin = 2, chains = 1, iter = 2000, warmup = iter/2,
                        seed = NA, M = NULL,
                        cmp_version = c("fast", "reference"),
                        cmp_tail_tol = 1e-6, tail_action = c("warn", "error"), ...) {

  if (identical(distri, "cmp")) {
    .check_counts(y)
  } else if (!is.numeric(y) || any(y != as.integer(y))) {
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


  stan_data <- list(y = c(y), T = length(y), ff = n_pred)
  if (distri == "cmp") {
    .check_scalar(n_pred, "n_pred", 0, .Machine$integer.max, integer = TRUE)
    cmp_version <- match.arg(cmp_version)
    tail_action <- match.arg(tail_action)
    if (is.null(M)) M <- max(300L, max(y))
    .check_scalar(M, "M", max(1L, max(y)), .Machine$integer.max - 2L,
                  integer = TRUE)
    .check_scalar(cmp_tail_tol, "cmp_tail_tol", 0, 1)
    if (cmp_tail_tol <= 0 || cmp_tail_tol >= 1) {
      stop("'cmp_tail_tol' must lie strictly between 0 and 1.")
    }
    stan_data$M <- as.integer(M)
  }

  stan_name <- switch(
    paste0(mod_type, "_", distri),
    "zi_poi" = "ZIP-INAR1.stan",
    "zi_nb"  = "ZINB-INAR1.stan",
    "zi_gp"  = "ZIGP-INAR1.stan",
    "h_poi"  = "HP-INAR1.stan",
    "h_nb"   = "HNB-INAR1.stan",
    "zi_cmp" = if (cmp_version == "fast") "ZICMP-INAR1-fast.stan" else
      "ZICMP-INAR1.stan",
    stop("Unsupported mod_type/distri combination: ", mod_type, "/", distri)
  )

  stan_file <- system.file("stan", stan_name, package = "ZIHINAR1")
  if (stan_file == "") {
    stop("Stan model file not found: ", stan_name, ". ",
         "In the package source, it must be located at inst/stan/", stan_name, ".")
  }

  stan_fit <- rstan::stan(file = stan_file,
                          data = stan_data,
                          thin = thin,
                          chains = chains,
                          iter = iter,
                          warmup = warmup,
                          seed = seed, ...)

  if (distri == "cmp") {
    .check_cmp_tail(stan_fit, cmp_tail_tol, tail_action)
    attr(stan_fit, "zihinar_info") <- list(
      distri = distri, mod_type = mod_type, y = as.integer(y),
      n_pred = as.integer(n_pred), M = as.integer(M),
      cmp_version = cmp_version, cmp_tail_tol = cmp_tail_tol)
  }

  return(stan_fit)
}
