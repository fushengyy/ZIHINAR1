install_if_missing("matrixStats")
library(matrixStats)

get_likelihood = function(y, alpha, rho, theta, mod_type, distri) {
  #theta = c(lambda, phi) for nb distribution
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
                   dpois(y[t], lambda = theta[1], log = TRUE)))
      }
      else{
        mu[t] = exp(dbinom(0, size = y[t-1], prob = alpha, log = TRUE)) *
          exp(dbinom(0, size = 1, prob = rho, log = TRUE) +
                dpois(y[t], lambda = theta[1], log = TRUE))
      }
      for (j in 1:pp) {
        if (y[t]==j) {
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha, log = TRUE)) *
            (exp(dbinom(1, size = 1, prob = rho, log = TRUE)) +
               exp(dbinom(0, size = 1, prob =  rho, log = TRUE) +
                     dpois(y[t]-j, lambda = theta[1], log = TRUE)))
        }
        else{
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha, log = TRUE)) *
            exp(dbinom(0, size = 1, prob = rho, log = TRUE) +
                  dpois(y[t]-j, lambda = theta[1], log = TRUE))
        }
      }
    }
  }

  if (mod_type == "zi" && distri == "nb") {
    betaBN = theta[2] / theta[1]
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
                   dnbinom(y[t], size = theta[2], prob = betaBN/(betaBN+1), log = TRUE)))
      }
      else{
        mu[t] = exp(dbinom(0, size = y[t-1], prob = alpha, log = TRUE)) *
          exp(dbinom(0, size = 1, prob = rho, log = TRUE) +
                dnbinom(y[t], size = theta[2], prob = betaBN/(betaBN+1), log = TRUE))
      }
      for (j in 1:pp) {
        if (y[t]==j) {
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha, log = TRUE)) *
            (exp(dbinom(1, size = 1, prob = rho, log = TRUE)) +
               exp(dbinom(0, size = 1, prob =  rho, log = TRUE) +
                     dnbinom(y[t]-j, size = theta[2], prob = betaBN/(betaBN+1), log = TRUE)))
        }
        else{
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha, log = TRUE)) *
            exp(dbinom(0, size = 1, prob = rho, log = TRUE) +
                  dnbinom(y[t]-j, size = theta[2], prob = betaBN/(betaBN+1), log = TRUE))
        }
      }
    }
  }

  if (mod_type == "h" && distri == "poi") {
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
          exp(log(1-rho) + dpois(y[t], lambda = theta[1], log = TRUE) -
                log(1-ppois(0, lambda = theta[1])))
      }
      for (j in 1:pp) {
        if (y[t]==j) {
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha, log = TRUE)) *
            exp(log(rho))
        }
        else{
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha, log = TRUE)) *
            exp(log(1-rho) + dpois(y[t]-j, lambda = theta[1], log = TRUE) -
                  log(1-ppois(0, lambda = theta[1])))
        }
      }
    }
  }

  if (mod_type == "h" && distri == "nb") {
    betaBN = theta[2] / theta[1]
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
          exp(log(1-rho) + dnbinom(y[t], size = theta[2], prob = betaBN/(betaBN+1), log = TRUE) -
                log(1-pnbinom(0, size = theta[2], prob = betaBN/(betaBN+1))))
      }
      for (j in 1:pp) {
        if (y[t]==j) {
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha, log = TRUE)) *
            exp(log(rho))
        }
        else{
          mu[t] = mu[t] + exp(dbinom(j, size = y[t-1], prob = alpha, log = TRUE)) *
            exp(log(1-rho) + dnbinom(y[t]-j, size = theta[2], prob = betaBN/(betaBN+1), log = TRUE) -
                  log(1-pnbinom(0, size = theta[2], prob = betaBN/(betaBN+1))))
        }
      }
    }
  }
  return(mu[2:length(y)])
}



get_loglik = function(y, alpha, rho, theta, mod_type, distri) {
  return(log(get_likelihood(y, alpha, rho, theta, mod_type, distri)))
}


mod_sel_criteria <- function(y, mod_type, distri, stanfit) {
  aic = data.frame(extract(stanfit, pars = "aic"))
  eaic = mean(aic[,1])

  bic = data.frame(extract(stanfit, pars = 'bic'))
  ebic = mean(bic[,1])

  if (distri == 'poi') {
    para.hat = summary(stanfit, pars = c("alpha", "rho", "lambda"))$summary
    alpha.hat = para.hat[1,1]
    rho.hat = para.hat[2,1]
    theta.hat = para.hat[3,1]
  }
  else {
    para.hat = summary(stanfit, pars = c("alpha", "rho", "lambda", "phi"))$summary
    alpha.hat = para.hat[1,1]
    rho.hat = para.hat[2,1]
    theta.hat = c(para.hat[3,1], para.hat[4,1])
  }

  logphat = sum(get_loglik(y, alpha = alpha.hat, rho = rho.hat, theta = theta.hat, mod_type, distri))
  ll = data.frame(extract(stanfit, pars = "ll"))
  pdic=2*(logphat - mean(ll[,1]))
  dic = -2*logphat + 2*pdic

  lik = data.frame(extract(stanfit, pars = "lik"))
  loglik = data.frame(extract(stanfit, pars = "log_lik"))
  pwaic1 = 2 * sum(log(colMeans(lik)) - colMeans(loglik))
  pwaic2 = sum(colVars(as.matrix(loglik)))
  lppd = sum(log(colMeans(lik)))
  waic1 = -2 * (lppd - pwaic1)
  waic2 = -2 * (lppd - pwaic2)

  all.criteria = data.frame(eaic, ebic, dic, waic1, waic2)
  colnames(all.criteria) = c("eaic", "ebic", "dic", "waic1", "waic2")
  return(all.criteria)
}
