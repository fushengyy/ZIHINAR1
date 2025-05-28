stanfit_est <- function(distri, stanfit) {
  if (distri == "poi") {
    qoi <- c("alpha", "rho", "lambda")
  } else {
    qoi <- c("alpha", "rho", "lambda", "phi")
  }

  est_all <- as.data.frame(summary(stanfit, pars = qoi)$summary)

  return(est_all)
}
