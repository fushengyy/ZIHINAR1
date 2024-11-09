/**
 * Time Series of Count Data: HP-INAR(1)
 * Log-Likelihood derived from:
 * ...
 */

data{
  int<lower=0> T;
  int y[T];
    int<lower=0> ff; // forecast
}
parameters{
  real<lower=0, upper=1> alpha;           // intercept
  real<lower=0> lambda;    // poisson parameter innovation
  real<lower=0, upper=1> rho;
}
transformed parameters{
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
      mu[t] = exp(binomial_lpmf(0| y[t-1], alpha))*exp(log(rho));
    }
    else{
      mu[t] = exp(binomial_lpmf(0| y[t-1], alpha))*exp(log1m(rho)+poisson_lpmf(y[t]|lambda)-poisson_lccdf(0|lambda));
    }
    for (j in 1:pp){
      if(y[t]==j){
        mu[t] += exp(binomial_lpmf(j| y[t-1], alpha))*exp(log(rho));
      }
      else{
        mu[t] += exp(binomial_lpmf(j| y[t-1], alpha))*exp(log1m(rho)+poisson_lpmf(y[t]-j|lambda)-poisson_lccdf(0|lambda));
      }
    }
  }
}

model{
  alpha  ~ uniform(0,1);
  rho ~   uniform(0,1);
  lambda  ~ lognormal(0, 2);
  for (t in 1:T) {
    target += log(mu[t]);
  }
}

generated quantities { // FOR PREDICTION
  // Generate posterior predictives
  int y_pred_ori[ff+1];
  int aa;
  int bb;
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
      bb=0;
      while(bb<1){
        bb = poisson_rng(lambda);
      }
      y_pred_ori[t] += bb;
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

  real aic = -2*ll+2*3;
  real bic = -2*ll+3*log(T-1);
}
