functions {
  real michaels_model(real r_a, real r_s, real r_l, real A_0, real beta, real rho, real time, real tau){
    return (time >= tau) * beta * (rho * ((exp(-r_s*(time-tau))-exp(-r_a*(time-tau))) / (r_a-r_s)) + (1-rho) * ((exp(-r_l*(time-tau))-exp(-r_a*(time-tau))) / (r_a-r_l)));
  }

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
  int<lower=0> Ny; // number of total data points
  int<lower=0> Nt; // number of total time points
  int<lower=0> K; // number of individuals
  vector<lower=0>[Ny] y; // antibody data
  vector<lower=0>[Nt] time; // time points
  int size_bt[K]; // number of exposures for each individual
  int size_y[K]; // number of antibody points for each individual 
  int size_T[K]; // number of time points for each individual
  vector<lower=0>[K] taus; // all exposure times
  real A_0;
  int Nyt;
  vector[Nyt] yrep_times;
}


parameters {
  // individual-level parameters
  vector<lower=0>[K] beta;
  vector<lower=0,upper=1>[K] rho;
  vector<lower=1>[K] d_a;
  vector<lower=1>[K] d_l;
  vector<lower=1>[K] d_s;
  
  real<lower=0> beta_sd;
  real<lower=0> rho_sd;
  real<lower=0> da_sd;
  real<lower=0> dl_sd;
  real<lower=0> ds_sd;
  
  //hyperparameters
  real<lower=0> beta_mu;
  real<lower=0,upper=1> rho_mu;
  real<lower=0> da_mu;
  real<lower=0> dl_mu;
  real<lower=0> ds_mu;

  real<lower=0> rate;
}


transformed parameters {
  vector<lower=0>[K] r_a;
  vector<lower=0>[K] r_l;
  vector<lower=0>[K] r_s;
  
  r_a = log(2)./d_a;
  r_s = log(2)./d_s;
  r_l = log(2)./d_l;
}


model {
  // individual-level parameters
  rho ~ lognormal(lognormal_mean(rho_mu, rho_sd), rho_sd);
  d_l ~ lognormal(lognormal_mean(dl_mu, dl_sd), dl_sd);
  d_a ~ lognormal(lognormal_mean(da_mu, da_sd), da_sd);
  d_s ~ lognormal(lognormal_mean(ds_mu, ds_sd), ds_sd);
  beta ~ lognormal(lognormal_mean(beta_mu, beta_sd), beta_sd);

  // hyperpriors
  rho_mu ~ beta(5, 2);
  dl_mu ~ lognormal(lognormal_mean(400, 0.2), 0.2);
  da_mu ~ lognormal(lognormal_mean(21, 0.1), 0.1);
  ds_mu ~ lognormal(lognormal_mean(3.5, 0.15), 0.15);
  beta_mu ~ lognormal(lognormal_mean(1500, 0.3), 0.3);
  
  rho_sd ~ lognormal(0, 1);
  dl_sd ~ lognormal(0, 1);
  da_sd ~ lognormal(0, 1);
  ds_sd ~ lognormal(0, 1);
  beta_sd ~ lognormal(0, 1);
  
  rate ~ lognormal(0, 1);
  
  for (i in 1:K) { // for each individual 
    vector[size_T[i]] time_seg = segment(time, sum(size_T[:i-1])+1, size_T[i]);  
    vector[size_T[i]] ab = A_0*exp(-r_l[i]*time_seg);
  
    for (j in 1:num_elements(time_seg)) { // loop through each time point
      ab[j] += michaels_model(r_a[i], r_s[i], r_l[i], A_0, beta[i], rho[i], time_seg[j], taus[i]);
    }
    segment(y, sum(size_T[:i-1])+1, size_T[i]) ~ gamma(vgamma_alpha(ab, rate), rate);
  }
}


generated quantities {
  int pos = 1;
  vector[Nyt*K] y_rep;
  vector[Nyt*K] y_rep_model;

  for (i in 1:K) {
    vector[Nyt] ab = A_0*exp(-r_l[i]*yrep_times);

    for (j in 1:Nyt) { // loop through each time point
      ab[j] += michaels_model(r_a[i], r_s[i], r_l[i], A_0, beta[i], rho[i], yrep_times[j], taus[i]);

      y_rep_model[pos] = ab[j];
      y_rep[pos] = gamma_rng(gamma_alpha(ab[j], rate), rate);
      pos += 1;
    }
  }
}

