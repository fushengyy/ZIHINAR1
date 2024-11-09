#' Simulate Time Series Data for ZI-INAR(1) or H-INAR(1) Models
#'
#' This function simulates time series data for Zero-Inflated INAR(1) (ZI-INAR(1)) or Hurdle INAR(1) (H-INAR(1)) models, using either Poisson or Negative Binomial distributions.
#'
#' @param n Integer specifying the number of observations to simulate.
#' @param alpha Numeric value between 0 and 1 representing the autoregressive parameter.
#' @param rho Numeric value between 0 and 1 representing the zero-inflation or hurdle parameter.
#' @param theta Numeric vector of model parameters. For Poisson, it should be \code{theta[1] = lambda}. For Negative Binomial, it should be \code{theta = c(lambda, phi)}, where \code{phi} is the size parameter.
#' @param mod_type Character string indicating the model type. Use \code{"zi"} for zero-inflated models and \code{"h"} for hurdle models.
#' @param distri Character string specifying the distribution. Options are \code{"poi"} for Poisson or \code{"nb"} for Negative Binomial.
#'
#' @return A numeric vector containing the simulated time series data.
#'
#' @references
#' Part of the implementation of this function was adapted from the \strong{ZINAR1} package.
#' The \strong{ZINAR1} package simulates first-order integer-valued autoregressive processes
#' with zero-inflated innovations (ZINAR(1)) and estimates its parameters under a frequentist approach.
#'
#' For more information about the ZINAR1 package, please refer to:
#'
#' Aldo M. Garay, João Vitor Ribeiro (2022).
#' \emph{ZINAR1: Simulates ZINAR(1) Model and Estimates Its Parameters Under Frequentist Approach}.
#' R package version 0.1.0. Available at: \url{https://CRAN.R-project.org/package=ZINAR1}.
#'
#' Garay, A. M., Ribeiro, J. V. (2021). First-Order Integer Valued AR Processes with Zero-Inflated Innovations.
#' In: \emph{Nonstationary Systems: Theory and Applications}, Springer.
#' DOI: \url{https://doi.org/10.1007/978-3-030-82110-4_2}.
#'
#' We acknowledge the original authors, Aldo M. Garay and João Vitor Ribeiro, for their contributions.
#'
#' @examples
#' \dontrun{
#'   # Simulate 100 observations from a Zero-Inflated Poisson INAR(1) model
#'   y_data <- data_simu(n = 100, alpha = 0.5, rho = 0.3, theta = c(5), mod_type = "zi", distri = "poi")
#'
#'   # Explore the simulated data
#'   if(!require('ZINARp')) {
#'     install.packages('ZINARp')
#'     library('ZINARp')
#'    }
#'   explore_zinarp(y_data)
#' }
#'
#' @importFrom VGAM rzipois rzinegbin
#' @importFrom actuar rztpois rztnbinom
#' @export
data_simu <- function(n, alpha, rho, theta, mod_type, distri) {
  # theta = c(lambda, phi) for nb distribution

  install_if_missing("VGAM")
  install_if_missing("actuar")

  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n)) {
    stop("The parameter 'n' must be a positive integer specifying the number of observations to simulate.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("The parameter 'alpha' must be a numeric value between 0 and 1 (excluding endpoints).")
  }
  if (!is.numeric(rho) || length(rho) != 1 || rho < 0 || rho > 1) {
    stop("The parameter 'rho' must be a numeric value between 0 and 1 (inclusive).")
  }
  if (!(mod_type %in% c("zi", "h"))) {
    stop("The parameter 'mod_type' must be either 'zi' for zero-inflated or 'h' for hurdle.")
  }
  if (!(distri %in% c("poi", "nb"))) {
    stop("The parameter 'distri' must be either 'poi' for Poisson or 'nb' for Negative Binomial.")
  }
  if (distri == "poi") {
    if (!is.numeric(theta) || length(theta) != 1 || theta <= 0) {
      stop("For Poisson distribution ('distri' = 'poi'), 'theta' must be a numeric vector of length 1 with a positive value.")
    }
  } else if (distri == "nb") {
    if (!is.numeric(theta) || length(theta) != 2) {
      stop("For Negative Binomial distribution ('distri' = 'nb'), 'theta' must be a numeric vector of length 2: c(lambda, phi).")
    }
    if (any(theta <= 0)) {
      stop("For Negative Binomial distribution ('distri' = 'nb'), both 'lambda' and 'phi' in 'theta' must be positive.")
    }
  }

  n_neg = 100

  y = s = vector()

  if (mod_type == "zi" && distri == "poi") {
    muy = theta[1] * (1 - rho) / (1 - alpha)
    y[1] = round(muy)
    w = VGAM::rzipois(n + n_neg, theta[1], rho)
  }

  if (mod_type == "zi" && distri == "nb") {
    muy = theta[1] * (1 - rho) / (1 - alpha)
    y[1] = round(muy)
    w = VGAM::rzinegbin(n + n_neg, pstr0 = rho, munb = theta[1], size = theta[2])
  }

  if (mod_type == "h" && distri == "poi") {
    muy = theta[1] * (1 - rho) / ((1 - alpha) * (1 - exp(-theta[1])))
    y[1] = round(muy)
    w <- runif(n + n_neg)
    w[w < rho] <- 0
    w[w > 0] <- actuar::rztpois(n + n_neg - length(w[w < rho]), theta[1])
  }

  if (mod_type == "h" && distri == "nb") {
    muy = theta[1] * (1 - rho) / ((1 - alpha) * (1 - (theta[2] / (theta[1] + theta[2]))^theta[2]))
    y[1] = round(muy)
    w <- runif(n + n_neg)
    w[w < rho] <- 0
    w[w > 0] <- actuar::rztnbinom(n + n_neg - length(w[w < rho]), theta[2], 1 - theta[1] / (theta[1] + theta[2]))
  }

  for (i in 2:(n + n_neg)) {
    s[i] = rbinom(1, y[i - 1], alpha)
    y[i] = s[i] + w[i]
  }

  y = y[-(1:n_neg)]

  return(y)
}
