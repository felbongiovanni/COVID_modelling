functions {
  vector vgamma_alpha(vector mean_ab, real beta) {
    return mean_ab*beta;
  }
  
  real gamma_alpha(real mean_ab, real beta) {
    return mean_ab*beta;
  }
  
  real lognormal_mean(real mu, real sigma) {
    return log(mu) - pow(sigma, 2)/2;
  }
}


data {
  int N;
  int K;
  int J;
  vector[N] ab;
  vector[N] times;
  int no_segments[K];
  int size_segments[J];
}


parameters {
  vector<lower=1>[K] d_lambda;
  
  real<lower=1> lambda_mu;
  real<lower=0> lambda_sd;

  real<lower=0> rate;
}

transformed parameters {
  vector[K] lambda = log(2)./d_lambda;
}

model {
  int pos = 1;
  d_lambda ~ lognormal(lognormal_mean(lambda_mu, lambda_sd), lambda_sd);
  
  lambda_sd ~ lognormal(0, 1);
  lambda_mu ~ lognormal(lognormal_mean(200, 0.3), 0.3);
  
  rate ~ lognormal(0, 1);
  
  for (i in 1:K) {
    for (j in 1:no_segments[i]) {
      vector[size_segments[pos]] time_seg = segment(times, sum(size_segments[:pos-1])+1, size_segments[pos]);
      vector[size_segments[pos]] ab_seg = segment(ab, sum(size_segments[:pos-1])+1, size_segments[pos]);
      
      vector[size_segments[pos]] ab_decay = ab_seg[1]*exp(-lambda[i]*time_seg);

      ab_seg ~ gamma(vgamma_alpha(ab_decay, rate), rate);
      pos += 1;
    }
  }
}


generated quantities {
  int pos = 1;
  int y_pos = 1;
  vector[N] y_rep;

  for (i in 1:K) {
    for (j in 1:no_segments[i]) {
      vector[size_segments[pos]] time_seg = segment(times, sum(size_segments[:pos-1])+1, size_segments[pos]);
      vector[size_segments[pos]] ab_seg = segment(ab, sum(size_segments[:pos-1])+1, size_segments[pos]);

      vector[size_segments[pos]] ab_decay = ab_seg[1]*exp(-/lambda[i]*time_seg);

      for (k in 1:size_segments[pos]) {
        y_rep[y_pos] = gamma_rng(gamma_alpha(ab_decay[k], rate), rate);
        y_pos += 1;
      }
      pos += 1;
    }
  }
}

