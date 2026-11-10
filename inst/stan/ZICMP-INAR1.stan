/**
 * Time Series of Count Data: ZICMP-INAR(1)
 * Binomial thinning and zero-inflated CMP innovations.
 * Conditional likelihood: y[2:T] given y[1].
 * CMP is approximated on 0:M in BOTH likelihood and forecast RNG.
 * Check cmp_tail_bound after fitting and refit with larger M if necessary.
 * Rate/dispersion priors match the supplied fast model (df = 5).
 */
 // Reference implementation: recompute CMP PMF for clarity.
functions {
  // Upper bound on omitted mass of the untruncated CMP distribution.
  // For k >= M+1, successive weights have ratio at most q=lambda/(M+2)^nu.
  // If q < 1, tail <= w[M+1]/(1-q). Convert this to a bound on tail / Z.
  real cmp_tail_bound_fn(real log_lambda, real nu, int M, real log_Z) {
    real log_q = log_lambda - nu * log(M + 2.0);
    if (log_q >= 0) return 1;
    {
      real log_tail = (M + 1.0) * log_lambda
                      - nu * lgamma(M + 2.0) - log1m_exp(log_q);
      return inv_logit(log_tail - log_Z);
    }
  }
  real cmp_log_Z(real lambda, real nu, int M) {
    vector[M + 1] log_weights;
    for (k in 0:M) {
      log_weights[k + 1] = k * log(lambda) - nu * lgamma(k + 1.0);
    }
    return log_sum_exp(log_weights);
  }

  real cmp_lpmf(int y, real lambda, real nu, int M) {
    if (y < 0 || y > M) return negative_infinity();
    return y * log(lambda) - nu * lgamma(y + 1.0)
           - cmp_log_Z(lambda, nu, M);
  }

  real zicmp_inar_lpmf(int y, int y_prev, real alpha,
                       real lambda, real nu, real rho, int M) {
    int p = min(y_prev, y);
    vector[p + 1] log_terms;
    for (j in 0:p) {
      int innovation = y - j;
      real log_innovation;
      if (innovation == 0) {
        log_innovation = log_sum_exp(log(rho), log1m(rho)
                                     + cmp_lpmf(0 | lambda, nu, M));
      } else {
        log_innovation = log1m(rho) + cmp_lpmf(innovation | lambda, nu, M);
      }
      log_terms[j + 1] = binomial_lpmf(j | y_prev, alpha) + log_innovation;
    }
    return log_sum_exp(log_terms);
  }

  int cmp_rng(real lambda, real nu, int M) {
    vector[M + 1] log_weights;
    for (k in 0:M) {
      log_weights[k + 1] = k * log(lambda) - nu * lgamma(k + 1.0);
    }
    return categorical_logit_rng(log_weights) - 1;
  }
}
data {
  int<lower=2> T;
  array[T] int<lower=0> y;
  int<lower=0> ff; // Number of FUTURE time points, excluding y[T].
  int<lower=1> M;  // CMP innovation cutoff; must also cover observed counts.
}
transformed data {
  if (M < max(y)) {
    reject("M must be at least max(y). Increase M; also check cmp_tail_bound.");
  }
}
parameters {
  real<lower=0, upper=1> alpha; // Binomial thinning probability.
  real<lower=0> lambda;        // CMP rate, not generally its mean.
  real<lower=0> nu;            // CMP dispersion; R data_simu uses disp.
  real<lower=0, upper=1> rho;   // Structural-zero probability of innovations.
}
transformed parameters {
  vector[T - 1] log_mu;
  for (t in 2:T) {
    log_mu[t - 1] = zicmp_inar_lpmf(y[t] | y[t - 1], alpha,
                                   lambda, nu, rho, M);
  }
}
model {
  alpha ~ uniform(0, 1);
  rho ~ uniform(0, 1);
  lambda ~ student_t(5, 0, 5); // Positive support gives a half-Student-t prior.
  nu ~ student_t(5, 0, 5);
  target += sum(log_mu);
}
generated quantities {
  // Package-compatible output shapes, shared by reference and fast models.
  array[ff] int y_pred;
  vector[T - 1] log_lik = log_mu;
  vector[T - 1] lik = exp(log_lik);
  real ll = sum(log_lik);
  // Compatibility scores evaluated at EACH draw, not classical MLE AIC/BIC.
  real aic = -2 * ll + 2 * 4;
  real bic = -2 * ll + 4 * log(T - 1.0);
  // Legacy field: sum of transition probabilities, NOT the joint likelihood.
  real likelihood = sum(lik);
  real<lower=0, upper=1> cmp_tail_bound;
  {
    real log_Z = cmp_log_Z(lambda, nu, M);
    cmp_tail_bound = cmp_tail_bound_fn(log(lambda), nu, M, log_Z);
  }
  if (ff > 0) {
    int previous_y = y[T];
    for (h in 1:ff) {
      int next_y = binomial_rng(previous_y, alpha);
      if (bernoulli_rng(rho) == 0) {
        next_y += cmp_rng(lambda, nu, M);
      }
      y_pred[h] = next_y;
      previous_y = next_y;
    }
  }
}
