/**
 * Time Series of Count Data: ZIGP-INAR(1)
 * Used the pmf given in the library HMMpa (dgenpois)
 * an also the paper by Harry Joe, Rong Zhu (2005) Bioemtrica Journal
 * https://onlinelibrary.wiley.com/doi/epdf/10.1002/bimj.200410102
 */

functions {
  real genpoisson_lpmf(int y, real lambda, real phi) {
    if (y < 0) return negative_infinity();
    real lam = lambda + phi * y;
    if (lam <= 0) return negative_infinity();
    return log(lambda) + (y - 1) * log(lam) - lam - lgamma(y + 1);
  }

  int generalized_poisson_rng(real lambda, real phi) {
    real u = uniform_rng(0, 1);
    real sum = 0;
    //int y;
    for (y in 0:999) {
      real p = exp(genpoisson_lpmf(y | lambda, phi));
      sum += p;
      if (u < sum) {
        return y;
      }
    }
    reject("RNG failed to converge within 1000 iterations.");
    return -1;  // this line is never reached, but required for return type
  }
}

data{
  int<lower=0> T;
  array[T] int y;
    int<lower=0> ff; // forecast
}
parameters{
  real<lower=0, upper=1> alpha;   // thinning parameter
  real<lower=0> lambda;    // GP parameter lambda
  real<lower=0,upper=1> phi;    //GP parameter eta

  real<lower=0, upper=1> rho; // probability of zero inflation
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
            mu[t] = exp(binomial_lpmf(0| y[t-1], alpha))*exp(log_sum_exp(bernoulli_lpmf(1 | rho), bernoulli_lpmf(0 | rho)
                      + genpoisson_lpmf(y[t]| lambda, phi)));
           }
           else{
            mu[t] = exp(binomial_lpmf(0| y[t-1], alpha))*exp(bernoulli_lpmf(0 | rho)+genpoisson_lpmf(y[t]| lambda, phi));
            }
        for (j in 1:pp){
             if(y[t]==j){
            mu[t] += exp(binomial_lpmf(j| y[t-1], alpha))*exp(log_sum_exp(bernoulli_lpmf(1 | rho), bernoulli_lpmf(0 | rho)
                      + genpoisson_lpmf(y[t]-j| lambda, phi)));
           }
           else{
            mu[t] += exp(binomial_lpmf(j| y[t-1], alpha))*exp(bernoulli_lpmf(0 | rho)+genpoisson_lpmf(y[t]-j| lambda, phi));
            }
        }
    }
}
model{
    alpha  ~ uniform(0,1);
    rho ~   uniform(0,1);
  lambda  ~ lognormal(0, 2);
  phi  ~ uniform(0,1);

for (t in 1:T) {
               target += log(mu[t]);
               }
}

generated quantities { // FOR PREDICTION
    // Generate posterior predictives
    array[ff+1] int y_pred_ori;
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
            y_pred_ori[t] +=generalized_poisson_rng(lambda, phi);
            }
    }

    array[ff] int y_pred;

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
