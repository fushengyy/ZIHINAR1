stanfit_est <- function(distri, stanfit) {
  if (distri == "poi") {
    qoi <- c("alpha", "rho", "lambda")
  } else if (distri == "cmp") {
    qoi <- c("alpha", "rho", "lambda", "nu")
  } else {
    qoi <- c("alpha", "rho", "lambda", "phi")
  }

  est_all <- as.data.frame(rstan::summary(stanfit, pars = qoi)$summary)

  return(est_all)
}
