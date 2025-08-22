lognormal_mean <- function(mu, sigma) {
  return(log(mu) - sigma^2 / 2)
}

lognormal_sd <- function(mu, sigma) {
  sqrt((exp(sigma^2) - 1) * exp(2 * mu + sigma^2))
}

antibody_model <- function(time, t_boost, beta, rho, d_a, d_l, d_s, A_0) {
  N <- length(t_boost)
  r_a <- log(2) / d_a
  r_l <- log(2) / d_l
  r_s <- log(2) / d_s

  ab <- A_0 * exp(-r_a * time)
  for (i in 1:N) {
    ab <- ab +
      as.numeric(time >= t_boost[i]) *
        beta[i] *
        (rho *
          ((exp(-r_s * (time - t_boost[i])) - exp(-r_a * (time - t_boost[i]))) /
            (r_a - r_s)) +
          (1 - rho) *
            ((exp(-r_l * (time - t_boost[i])) -
              exp(-r_a * (time - t_boost[i]))) /
              (r_a - r_l)))
  }
  return(ab)
}

michaels_antibody_model <- function(
  time,
  t_boost,
  beta,
  rho,
  d_a,
  d_l,
  d_s,
  A_0
) {
  N <- length(t_boost)
  r_a <- log(2) / d_a
  r_l <- log(2) / d_l
  r_s <- log(2) / d_s

  ab <- A_0 * exp(-r_l * time)
  for (i in 1:N) {
    ab <- ab +
      as.numeric(time >= t_boost[i]) *
        beta[i] *
        (rho *
          ((exp(-r_s * (time - t_boost[i])) - exp(-r_a * (time - t_boost[i]))) /
            (r_a - r_s)) +
          (1 - rho) *
            ((exp(-r_l * (time - t_boost[i])) -
              exp(-r_a * (time - t_boost[i]))) /
              (r_a - r_l)))
  }
  return(ab)
}

michaels_antibody_model_pp <- function(
  time,
  t_boost,
  beta,
  rho,
  d_a,
  d_l,
  d_s,
  A_0
) {
  N <- length(t_boost)
  r_a <- log(2) / d_a
  r_l <- log(2) / d_l
  r_s <- log(2) / d_s

  ab <- A_0 * exp(-r_l * time)
  for (i in 1:N) {
    ab <- ab +
      as.numeric(time >= t_boost[i]) *
        beta *
        (rho *
          ((exp(-r_s * (time - t_boost[i])) - exp(-r_a * (time - t_boost[i]))) /
            (r_a - r_s)) +
          (1 - rho) *
            ((exp(-r_l * (time - t_boost[i])) -
              exp(-r_a * (time - t_boost[i]))) /
              (r_a - r_l)))
  }
  return(ab)
}

michaels_antibody_model_vector <- function(
  time,
  t_boost,
  beta,
  rho,
  d_a,
  d_l,
  d_s,
  A_0,
  rate
) {
  N <- length(t_boost)
  r_a <- log(2) / d_a
  r_l <- log(2) / d_l
  r_s <- log(2) / d_s

  ab <- A_0 * exp(-r_l * time)
  for (i in 1:N) {
    ab <- ab +
      as.numeric(time >= t_boost[i]) *
        beta *
        (rho *
          ((exp(-r_s * (time - t_boost[i])) - exp(-r_a * (time - t_boost[i]))) /
            (r_a - r_s)) +
          (1 - rho) *
            ((exp(-r_l * (time - t_boost[i])) -
              exp(-r_a * (time - t_boost[i]))) /
              (r_a - r_l)))
  }
  ab_noise <- rgamma(length(ab), shape = ab * rate, rate = rate)
  return(ab_noise)
}

michaels_antibody_model_noise <- function(
  time,
  t_boost,
  beta,
  rho,
  d_a,
  d_l,
  d_s,
  A_0,
  rate
) {
  N <- length(t_boost)
  r_a <- log(2) / d_a
  r_l <- log(2) / d_l
  r_s <- log(2) / d_s

  ab <- A_0 * exp(-r_l * time)
  for (i in 1:N) {
    ab <- ab +
      as.numeric(time >= t_boost[i]) *
        beta[i] *
        (rho *
          ((exp(-r_s * (time - t_boost[i])) - exp(-r_a * (time - t_boost[i]))) /
            (r_a - r_s)) +
          (1 - rho) *
            ((exp(-r_l * (time - t_boost[i])) -
              exp(-r_a * (time - t_boost[i]))) /
              (r_a - r_l)))
  }
  ab_noise <- rgamma(length(ab), shape = ab * rate, rate = rate)
  return(ab_noise)
}

generating_data <- function(
  A0_mu,
  beta_mu,
  rho_mu,
  da_mu,
  ds_mu,
  dl_mu,
  num,
  t_boost,
  rate
) {
  stan_df <- list()
  priors_df <- list()
  for (i in 1:num) {
    A0_prior <- 1
    beta_prior <- rlnorm(1, lognormal_mean(beta_mu, 0.05), 0.05)
    rho_prior <- rlnorm(1, lognormal_mean(rho_mu, 0.05), 0.05)
    ds_prior <- rlnorm(1, lognormal_mean(ds_mu, 0.05), 0.05)
    dl_prior <- rlnorm(1, lognormal_mean(dl_mu, 0.07), 0.07)
    da_prior <- rlnorm(1, lognormal_mean(da_mu, 0.05), 0.05)
    prim_inf <- sample(-150:-10, 1)
    if (abs(prim_inf) < 25) {
      times <- c(
        0,
        floor(seq(abs(prim_inf), sample(500:600, 1), length.out = 10))
      )
      dd <- times[3]
    } else {
      times <- c(
        0,
        21,
        floor(seq(abs(prim_inf), sample(500:600, 1), length.out = 9))
      )
      dd <- prim_inf
    }
    A <- michaels_antibody_model(
      times,
      t_boost,
      beta_prior,
      rho_prior,
      da_prior,
      dl_prior,
      ds_prior,
      A0_prior
    )
    A <- rgamma(length(A), shape = A * rate, rate = rate)

    stan_df[[i]] <- list(
      A,
      A[2],
      A[3:length(A)],
      times,
      times[3:length(times)],
      i,
      A[3]
    )
    priors_df[[i]] <- list(
      A0_prior,
      t_boost,
      beta_prior,
      rho_prior,
      ds_prior,
      dl_prior,
      da_prior,
      rate,
      i,
      dd,
      length(t_boost),
      length(A),
      length(times)
    )
  }
  return(list(stan_df, priors_df))
}


generating_data_multi_boost <- function(
  A0_mu,
  beta_mu,
  rho_mu,
  da_mu,
  ds_mu,
  dl_mu,
  num,
  taus,
  times,
  rate
) {
  stan_df <- list()
  priors_df <- list()
  for (i in 1:num) {
    rho_prior <- rlnorm(1, lognormal_mean(rho_mu, 0.05), 0.05)
    ds_prior <- rlnorm(1, lognormal_mean(ds_mu, 0.2), 0.2)
    dl_prior <- rlnorm(1, lognormal_mean(dl_mu, 0.3), 0.3)
    da_prior <- rlnorm(1, lognormal_mean(da_mu, 0.1), 0.1)

    ind_times <- times[[i]]
    ind_tau <- taus[[i]]

    beta_prior <- rlnorm(length(ind_tau), lognormal_mean(beta_mu, 0.3), 0.3)

    A <- michaels_antibody_model(
      ind_times,
      ind_tau,
      beta_prior,
      rho_prior,
      da_prior,
      dl_prior,
      ds_prior,
      A0_mu
    )
    A <- rgamma(length(A), shape = A * rate, rate = rate)

    stan_df[[i]] <- list(A, ind_times, i, ind_tau)
    priors_df[[i]] <- list(
      A0_mu,
      ind_tau,
      beta_prior,
      rho_prior,
      ds_prior,
      dl_prior,
      da_prior,
      i,
      length(ind_tau),
      length(A),
      length(ind_times)
    )
  }
  return(list(stan_df, priors_df))
}


bounded_rlnorm <- function(n, mean, sd, lower, upper) {
  samples <- rlnorm(n, lognormal_mean(mean, sd), sd)
  idx <- which(samples < lower | samples > upper)

  while (length(idx) > 0) {
    samples[idx] <- rlnorm(length(idx), lognormal_mean(mean, sd), sd)
    idx <- which(samples < lower | samples > upper)
  }
  return(samples)
}
