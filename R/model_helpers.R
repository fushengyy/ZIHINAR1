# Internal helpers for CMP fitting and model selection.
.check_counts <- function(y) {
  if (!is.numeric(y) || is.complex(y) || !is.null(dim(y)) || length(y) < 2L ||
      any(!is.finite(y)) || any(y < 0 | y != floor(y)) ||
      any(y > .Machine$integer.max)) {
    stop("The parameter 'y' must be a vector of at least two finite, nonnegative integers.")
  }
  invisible(NULL)
}

.check_scalar <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x))) {
    stop("The parameter '", name, "' must be a finite ",
         if (integer) "integer" else "number", " in [", lower, ", ", upper, "].")
  }
  invisible(NULL)
}

.check_fit <- function(stan_fit, distri = NULL, mod_type = NULL, y = NULL) {
  if (!inherits(stan_fit, "stanfit") || stan_fit@mode != 0L) {
    stop("Provide a 'stanfit' object with retained posterior draws.")
  }
  info <- attr(stan_fit, "zihinar_info", exact = TRUE)
  if (!is.null(info)) {
    if (!is.null(distri) && !identical(distri, info$distri)) {
      stop("The requested distribution does not match the fitted model.")
    }
    if (!is.null(mod_type) && !identical(mod_type, info$mod_type)) {
      stop("The requested model type does not match the fitted model.")
    }
    if (!is.null(y) && !identical(as.integer(y), info$y)) {
      stop("The supplied data do not match the data used to fit the model.")
    }
  }
  invisible(info)
}

# Always return draws x components, including a single draw or component.
.draw_matrix <- function(stan_fit, variable) {
  x <- rstan::extract(stan_fit, pars = variable, permuted = TRUE,
                      inc_warmup = FALSE)[[variable]]
  if (is.null(x) || length(x) == 0L) stop("No posterior draws found for '", variable, "'.")
  if (is.null(dim(x))) return(matrix(x, ncol = 1L))
  matrix(x, nrow = dim(x)[1L])
}

.check_cmp_tail <- function(stan_fit, tol = 1e-6, action = c("warn", "error")) {
  action <- match.arg(action)
  bounds <- .draw_matrix(stan_fit, "cmp_tail_bound")
  if (any(!is.finite(bounds)) || any(bounds > tol)) {
    msg <- paste0("CMP truncation is inadequate (maximum tail bound = ",
                  signif(max(bounds), 4), "). Increase M and refit before using the results.")
    if (action == "error") stop(msg) else warning(msg, call. = FALSE)
  }
  invisible(max(bounds))
}

.cmp_log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

.cmp_loglik <- function(y, alpha, rho, lambda, nu, M) {
  .check_counts(y)
  .check_scalar(alpha, "alpha", 0, 1)
  .check_scalar(rho, "rho", 0, 1)
  .check_scalar(lambda, "lambda", 0)
  .check_scalar(nu, "nu", 0)
  if (lambda <= 0 || nu <= 0) stop("'lambda' and 'nu' must be positive.")
  if (is.null(M)) M <- max(300L, max(y))
  .check_scalar(M, "M", max(1L, max(y)), .Machine$integer.max - 2L, integer = TRUE)
  k <- 0:M
  log_weights <- k * log(lambda) - nu * lgamma(k + 1)
  log_innovation <- log1p(-rho) + log_weights - .cmp_log_sum_exp(log_weights)
  log_innovation[1L] <- .cmp_log_sum_exp(c(log(rho), log_innovation[1L]))
  vapply(seq.int(2L, length(y)), function(t) {
    j <- 0:min(y[t - 1L], y[t])
    .cmp_log_sum_exp(stats::dbinom(j, y[t - 1L], alpha, log = TRUE) +
                       log_innovation[y[t] - j + 1L])
  }, numeric(1))
}

.cmp_model_selection <- function(y, mod_type, stanfit, M) {
  if (mod_type != "zi") stop("CMP is supported only for zero-inflated models.")
  .check_counts(y)
  info <- .check_fit(stanfit, "cmp", mod_type, y)
  if (!is.null(M)) {
    .check_scalar(M, "M", max(1L, max(y)), .Machine$integer.max - 2L, integer = TRUE)
  }
  if (!is.null(info$M)) {
    if (!is.null(M) && M != info$M) stop("M must equal the cutoff used for fitting.")
    M <- info$M
  }
  if (is.null(M)) stop("Supply the exact CMP cutoff M for fits made directly with Stan.")
  .check_cmp_tail(stanfit, if (is.null(info$cmp_tail_tol)) 1e-6 else info$cmp_tail_tol)
  loglik <- .draw_matrix(stanfit, "log_lik")
  if (ncol(loglik) != length(y) - 1L || nrow(loglik) < 2L || any(!is.finite(loglik))) {
    stop("Expected at least two draws with T-1 finite transition log likelihoods.")
  }
  est <- stanfit_est("cmp", stanfit)
  logphat <- sum(.cmp_loglik(y, est["alpha", "mean"], est["rho", "mean"],
                           est["lambda", "mean"], est["nu", "mean"], M))
  Dbar <- -2 * mean(rowSums(loglik))
  log_mean_lik <- apply(loglik, 2L, .cmp_log_sum_exp) - log(nrow(loglik))
  lppd <- sum(log_mean_lik)
  pwaic1 <- 2 * sum(log_mean_lik - colMeans(loglik))
  pwaic2 <- sum(matrixStats::colVars(loglik))
  data.frame(eaic = Dbar + 2 * 4, ebic = Dbar + 4 * log(ncol(loglik)),
             dic = 2 * Dbar + 2 * logphat,
             waic1 = -2 * (lppd - pwaic1), waic2 = -2 * (lppd - pwaic2))
}
