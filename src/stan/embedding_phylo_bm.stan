// embedding_phylo_bm.stan
// Bayesian phylogenetic model for protein embedding evolution
//
// Brownian Motion of embedding dimensions along a phylogeny.
// For each embedding dimension k = 1..D:
//   z[, k] ~ MVN(z_anc[k] * 1_S, rate[k] * C)
//
// Non-centered parameterization:
//   z[, k] = z_anc[k] + sqrt(rate[k]) * L_C * z_raw[, k]
//
// The model jointly estimates:
//   - z_anc: the ancestral embedding at the root
//   - rate: per-dimension BM evolutionary rates
//   - z: per-species latent embeddings
//
// The observed embeddings have Gaussian measurement noise.

data {
  int<lower=1> S;                    // number of species (tips)
  int<lower=1> D;                    // embedding dimensionality

  matrix[S, D] z_obs;               // observed embedding matrix
  matrix[S, S] C;                   // phylogenetic VCV (from ape::vcv)

  real<lower=0> obs_sigma;          // measurement noise sd

  vector[D] z_anc_prior_mean;       // prior mean for ancestral embedding
  vector<lower=0>[D] z_anc_prior_sd; // prior sd for ancestral embedding
  real<lower=0> rate_scale;          // half-normal scale for BM rates
}

transformed data {
  cholesky_factor_cov[S] L_C = cholesky_decompose(C);
}

parameters {
  vector[D] z_anc;
  vector<lower=0>[D] rate;
  matrix[S, D] z_raw;
}

transformed parameters {
  matrix[S, D] z;
  for (k in 1:D) {
    z[, k] = z_anc[k] + sqrt(rate[k]) * (L_C * z_raw[, k]);
  }
}

model {
  z_anc ~ normal(z_anc_prior_mean, z_anc_prior_sd);
  rate ~ normal(0, rate_scale);
  for (k in 1:D) {
    z_raw[, k] ~ std_normal();
  }
  for (s in 1:S) {
    z_obs[s] ~ normal(z[s], obs_sigma);
  }
}

generated quantities {
  matrix[S, D] z_pred;
  for (s in 1:S) {
    for (k in 1:D) {
      z_pred[s, k] = normal_rng(z[s, k], obs_sigma);
    }
  }
}
