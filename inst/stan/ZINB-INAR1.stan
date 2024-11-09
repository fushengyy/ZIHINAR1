/**
 * Time Series of Count Data: ZINB-INAR(1)
 * Log-Likelihood derived from:
 * ??
 */

data{
  int<lower=0> T;
  int y[T];
  int<lower=0> ff; // forecast
}

parameters{
  real<lower=0, upper=1> alpha;       // intercept
  real<lower=0> lambda;    // negative Binomial parameter innovation: mean
  real<lower=0> phi;    // negative Binomial parameter innovation: dispersion
  real<lower=0, upper=1> rho; // zero inflated parameter
}
transformed parameters{
  real<lower=0> betaBN;
  betaBN=phi/lambda;
  vector[T] mu;
  if (y[1]==0){
    mu[1]=1;
  }
  else{
    mu[1] = y[1];
  }
  for (t in 2:T){
    int pp=min(y[t-1:t]);
    if(y[t]==0){
      mu[t] = exp(binomial_lpmf(0| y[t-1], alpha))*exp(log_sum_exp(bernoulli_lpmf(1 | rho), bernoulli_lpmf(0 | rho)
                      + neg_binomial_lpmf(y[t] | phi, betaBN)));
    }
    else{
      mu[t] = exp(binomial_lpmf(0| y[t-1], alpha))*exp(bernoulli_lpmf(0 | rho)+neg_binomial_lpmf(y[t]|phi, betaBN));
    }
    for (j in 1:pp){
      if(y[t]==j){
        mu[t] += exp(binomial_lpmf(j| y[t-1], alpha))*exp(log_sum_exp(bernoulli_lpmf(1 | rho), bernoulli_lpmf(0 | rho)
                      + neg_binomial_lpmf(y[t]-j | phi, betaBN)));
      }
      else{
        mu[t] += exp(binomial_lpmf(j| y[t-1], alpha))*exp(bernoulli_lpmf(0 | rho)+neg_binomial_lpmf(y[t]-j|phi, betaBN));
      }
    }
  }
}

model{
  alpha  ~ uniform(0,1);
  rho ~   uniform(0,1);
  lambda  ~ lognormal(0, 2);
  phi  ~ lognormal(0, sqrt(2));

  for (t in 1:T) {
    target += log(mu[t]);
  }
}

generated quantities { // FOR PREDICTION
  // Generate posterior predictives
  int y_pred_ori[ff+1];
  int aa;

  // First P points are known
  y_pred_ori[1] = y[T];

  // Posterior predictive
  for (t in 2:ff+1){
    y_pred_ori[t] = binomial_rng(y_pred_ori[t-1], alpha);
    aa = bernoulli_rng(rho);
    if(aa==1){
      y_pred_ori[t] += 0;
    }
    else{
      y_pred_ori[t] +=neg_binomial_rng(phi, betaBN);
    }
  }

  int y_pred[ff];

  for (t in 1:ff) {
    y_pred[t] = y_pred_ori[t+1];
  }

  vector[T-1] log_lik;
  vector[T-1] lik;

  for (t in 2:T) {
    log_lik[t-1] = log(mu[t]);
    lik[t-1] = mu[t];
  }

  real ll = sum(log_lik);
  real likelihood = sum(lik);

  real aic = -2*ll+2*4;
  real bic = -2*ll+4*log(T-1);
}

