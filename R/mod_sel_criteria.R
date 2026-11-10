#' @importFrom stats dbinom dnbinom dpois median pnbinom ppois
#' @import matrixStats
#' @import HMMpa

dgenpois_log <- function(y, lambda, phi) {
  y <- as.integer(y)
  out <- rep(-Inf, length(y))

  ok_y <- !is.na(y) & y >= 0L
  if (any(ok_y)) {
    yy      <- y[ok_y]
    lam_y   <- lambda + phi * yy
    ok_lam  <- is.finite(lam_y) & (lam_y > 0)
    if (any(ok_lam)) {
      idx    <- which(ok_y)[ok_lam]
      yy_ok  <- yy[ok_lam]
      lam_ok <- lam_y[ok_lam]
      out[idx] <- log(lambda) + (yy_ok - 1) * log(lam_ok) - lam_ok - lgamma(yy_ok + 1)
    }
  }
  out
}

get_likelihood = function(y, alpha, rho, lambda, disp = NA_real_, mod_type, distri, M = NULL) {

  if (distri == "cmp") {
    return(exp(get_loglik(y, alpha, rho, lambda, disp, mod_type, distri, M)))
  }

  if (mod_type == "zi" && distri == "poi") {
    mu = c()
    if (y[1]==0){
      mu[1] = 1
    }
    else{
      mu[1] = y[1]
    }
    for (t in 2:length(y)) {
      pp = min(y[(t-1):t])
      if (y[t]==0) {
        mu[t] = exp(dbinom(0, size = y[t-1], prob = alpha, log = TRUE)) *
          (exp(dbinom(1, size = 1, prob = rho, log = TRUE)) +
             exp(dbinom(0, size = 1, prob =  rho, log = TRUE) +
                   dpois(y[t], lambda = lambda, log = TRUE)))
      }
      else{
        mu[t] = exp(dbinom(0, size = y[t-1], prob = alpha, log = TRUE)) *
          exp(dbinom(0, size = 1, prob = rho, log = TRUE) +
                dpois(y[t], lambda = lambda, log = TRUE))
      }
      for (j in seq_len(pp)) {
        if (y[t]==j) {
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha,
                                     log = TRUE)) *
            (exp(dbinom(1, size = 1, prob = rho, log = TRUE)) +
               exp(dbinom(0, size = 1, prob =  rho, log = TRUE) +
                     dpois(y[t]-j, lambda = lambda, log = TRUE)))
        }
        else{
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha,
                                     log = TRUE)) *
            exp(dbinom(0, size = 1, prob = rho, log = TRUE) +
                  dpois(y[t]-j, lambda = lambda, log = TRUE))
        }
      }
    }
  }

  else if (mod_type == "zi" && distri == "nb") {
    betaBN = disp / lambda
    mu = c()
    if (y[1]==0){
      mu[1] = 1
    }
    else{
      mu[1] = y[1]
    }
    for (t in 2:length(y)) {
      pp = min(y[(t-1):t])
      if (y[t]==0) {
        mu[t] = exp(dbinom(0, size = y[t-1], prob = alpha, log = TRUE)) *
          (exp(dbinom(1, size = 1, prob = rho, log = TRUE)) +
             exp(dbinom(0, size = 1, prob =  rho, log = TRUE) +
                   dnbinom(y[t], size = disp, prob = betaBN/(betaBN+1),
                           log = TRUE)))
      }
      else{
        mu[t] = exp(dbinom(0, size = y[t-1], prob = alpha, log = TRUE)) *
          exp(dbinom(0, size = 1, prob = rho, log = TRUE) +
                dnbinom(y[t], size = disp, prob = betaBN/(betaBN+1),
                        log = TRUE))
      }
      for (j in seq_len(pp)) {
        if (y[t]==j) {
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha,
                                     log = TRUE)) *
            (exp(dbinom(1, size = 1, prob = rho, log = TRUE)) +
               exp(dbinom(0, size = 1, prob =  rho, log = TRUE) +
                     dnbinom(y[t]-j, size = disp, prob = betaBN/(betaBN+1),
                             log = TRUE)))
        }
        else{
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha,
                                     log = TRUE)) *
            exp(dbinom(0, size = 1, prob = rho, log = TRUE) +
                  dnbinom(y[t]-j, size = disp, prob = betaBN/(betaBN+1),
                          log = TRUE))
        }
      }
    }
  }

  else if (mod_type == "h" && distri == "poi") {
    mu = c()
    if (y[1]==0){
      mu[1] = 1
    }
    else{
      mu[1] = y[1]
    }
    for (t in 2:length(y)) {
      pp = min(y[(t-1):t])
      if (y[t]==0) {
        mu[t] = exp(dbinom(0, size = y[t-1], prob = alpha, log = TRUE)) *
          exp(log(rho))
      }
      else{
        mu[t] = exp(dbinom(0, size = y[t-1], prob = alpha, log = TRUE)) *
          exp(log(1-rho) + dpois(y[t], lambda = lambda, log = TRUE) -
                log(1-ppois(0, lambda = lambda)))
      }
      for (j in seq_len(pp)) {
        if (y[t]==j) {
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha,
                                     log = TRUE)) * exp(log(rho))
        }
        else{
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha,
                                     log = TRUE)) *
            exp(log(1-rho) + dpois(y[t]-j, lambda = lambda, log = TRUE) -
                  log(1-ppois(0, lambda = lambda)))
        }
      }
    }
  }

  else if (mod_type == "h" && distri == "nb") {
    betaBN = disp / lambda
    mu = c()
    if (y[1]==0){
      mu[1] = 1
    }
    else{
      mu[1] = y[1]
    }
    for (t in 2:length(y)) {
      pp = min(y[(t-1):t])
      if (y[t]==0) {
        mu[t] = exp(dbinom(0, size = y[t-1], prob = alpha, log = TRUE)) *
          exp(log(rho))
      }
      else{
        mu[t] = exp(dbinom(0, size = y[t-1], prob = alpha, log = TRUE)) *
          exp(log(1-rho) + dnbinom(y[t], size = disp,
                                   prob = betaBN/(betaBN+1), log = TRUE) -
                log(1-pnbinom(0, size = disp, prob = betaBN/(betaBN+1))))
      }
      for (j in seq_len(pp)) {
        if (y[t]==j) {
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha,
                                     log = TRUE)) * exp(log(rho))
        }
        else{
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha,
                                     log = TRUE)) *
            exp(log(1-rho) + dnbinom(y[t]-j, size = disp,
                                     prob = betaBN/(betaBN+1), log = TRUE) -
                  log(1-pnbinom(0, size = disp, prob = betaBN/(betaBN+1))))
        }
      }
    }
  }

  else if (mod_type == "zi" && distri == "gp") {
    mu = c()
    if (y[1]==0){
      mu[1] = 1
    }
    else{
      mu[1] = y[1]
    }
    for (t in 2:length(y)) {
      pp = min(y[(t-1):t])
      if (y[t]==0) {
        mu[t] = exp(dbinom(0, size = y[t-1], prob = alpha, log = TRUE)) *
          (exp(dbinom(1, size = 1, prob = rho, log = TRUE)) +
             exp(dbinom(0, size = 1, prob =  rho, log = TRUE) +
                   dgenpois_log(y[t], lambda = lambda, phi = disp)))
      }
      else{
        mu[t] = exp(dbinom(0, size = y[t-1], prob = alpha, log = TRUE)) *
          exp(dbinom(0, size = 1, prob = rho, log = TRUE) +
                dgenpois_log(y[t], lambda = lambda, phi = disp))
      }
      for (j in seq_len(pp)) {
        if (y[t]==j) {
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha, log = TRUE)) *
            (exp(dbinom(1, size = 1, prob = rho, log = TRUE)) +
               exp(dbinom(0, size = 1, prob =  rho, log = TRUE) +
                     dgenpois_log(y[t]-j, lambda = lambda, phi = disp)))
        }
        else{
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha, log = TRUE)) *
            exp(dbinom(0, size = 1, prob = rho, log = TRUE) +
                  dgenpois_log(y[t]-j, lambda = lambda, phi = disp))
        }
      }
    }
  }

  else {
    stop("Unsupported mod_type/distri combination: ", mod_type, "/", distri)
  }

  return(mu[2:length(y)])
}



get_loglik <- function(y, alpha, rho, lambda, disp = NA_real_, mod_type, distri, M = NULL) {
  if (distri == "cmp") {
    if (mod_type != "zi") stop("CMP is supported only for zero-inflated models.")
    return(.cmp_loglik(y, alpha, rho, lambda, disp, M))
  }
  log(get_likelihood(y, alpha, rho, lambda, disp, mod_type, distri))
}


mod_sel_criteria <- function(y, mod_type, distri, stanfit, M = NULL) {
  if (distri == "cmp") {
    return(.cmp_model_selection(y, mod_type, stanfit, M))
  }

  aic = data.frame(extract(stanfit, pars = "aic"))
  eaic = mean(aic[,1])

  bic = data.frame(extract(stanfit, pars = 'bic'))
  ebic = mean(bic[,1])

  if (distri == 'poi') {
    para.hat = summary(stanfit, pars = c("alpha", "rho", "lambda"))$summary
    alpha.hat = para.hat[1,1]
    rho.hat = para.hat[2,1]
    lambda.hat = para.hat[3,1]

    logphat = sum(get_loglik(y, alpha = alpha.hat, rho = rho.hat,
                             lambda = lambda.hat, disp = NA,
                             mod_type, distri))
  }
  else {
    para.hat = summary(stanfit,
                       pars = c("alpha", "rho", "lambda", "phi"))$summary
    alpha.hat = para.hat[1,1]
    rho.hat = para.hat[2,1]
    lambda.hat = para.hat[3,1]
    disp.hat = para.hat[4,1]

    logphat = sum(get_loglik(y, alpha = alpha.hat, rho = rho.hat,
                             lambda = lambda.hat, disp = disp.hat,
                             mod_type, distri))
  }

  ll = data.frame(extract(stanfit, pars = "ll"))
  Dbar = -2 * mean(ll[,1])
  Dhat = -2 * logphat
  dic = 2 * Dbar - Dhat

  lik = data.frame(extract(stanfit, pars = "lik"))
  loglik = data.frame(extract(stanfit, pars = "log_lik"))
  pwaic1 = 2 * sum(log(colMeans(lik)) - colMeans(loglik))
  pwaic2 = sum(matrixStats::colVars(as.matrix(loglik)))
  lppd = sum(log(colMeans(lik)))
  waic1 = -2 * (lppd - pwaic1)
  waic2 = -2 * (lppd - pwaic2)

  all.criteria = data.frame(eaic, ebic, dic, waic1, waic2)
  colnames(all.criteria) = c("eaic", "ebic", "dic", "waic1", "waic2")
  return(all.criteria)
}
