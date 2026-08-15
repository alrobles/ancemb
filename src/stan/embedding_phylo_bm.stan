// embedding_phylo_bm.stan
// Bayesian phylogenetic model for protein embedding evolution
//
// Models embedding dimensions as continuous traits evolving under
// Brownian Motion along a phylogeny.
//
// For each embedding dimension k = 1..D:
//   z[, k] ~ MVN(z_anc[k] * 1_S, rate[k] * C)
//
// Non-centered parameterization for efficient HMC:
//   z[, k] = z_anc[k] + sqrt(rate[k]) * L_C * z_raw[, k]
//
// The observed embeddings are assumed to have Gaussian measurement noise,
// which regularizes the inference and accounts for pLM stochasticity.

data {
  int<lower=1> S;                    // number of species (tips)
  int<lower=1> D;                    // embedding dimensionality (after PCA)

  // Observed embeddings (one per species, after mean-pooling and PCA)
  matrix[S, D] z_obs;               // observed embedding matrix

  // Phylogenetic variance-covariance matrix (from ape::vcv)
  matrix[S, S] C;                   // phylogenetic VCV

  // Observation noise (estimated or fixed)
  real<lower=0> obs_sigma;          // measurement noise sd per dimension

  // Prior hyperparameters
  vector[D] z_anc_prior_mean;       // prior mean for ancestral embedding
  vector<lower=0>[D] z_anc_prior_sd; // prior sd for ancestral embedding
  real<lower=0> rate_scale;          // half-normal scale for BM rates
}

transformed data {
  // Pre-compute Cholesky of phylogenetic VCV (done once)
  cholesky_factor_cov[S] L_C = cholesky_decompose(C);
}

parameters {
  // Ancestral embedding (root state)
  vector[D] z_anc;

  // BM evolutionary rates (one per embedding dimension)
  vector<lower=0>[D] rate;

  // Non-centered species deviations (iid standard normal)
  matrix[S, D] z_raw;
}

transformed parameters {
  // Species embeddings via non-centered BM parameterization
  matrix[S, D] z;
  for (k in 1:D) {
    z[, k] = z_anc[k] + sqrt(rate[k]) * (L_C * z_raw[, k]);
  }
}

model {
  // --- Priors on ancestral embedding ---
  z_anc ~ normal(z_anc_prior_mean, z_anc_prior_sd);

  // --- Prior on BM rates (half-normal) ---
  rate ~ normal(0, rate_scale);

  // --- Non-centered BM deviations: iid N(0,1) ---
  for (k in 1:D) {
    z_raw[, k] ~ std_normal();
  }

  // --- Observation model: embeddings observed with Gaussian noise ---
  for (s in 1:S) {
    z_obs[s] ~ normal(z[s], obs_sigma);
  }
}

generated quantities {
  // Posterior predictive check: simulated observations
  matrix[S, D] z_pred;
  for (s in 1:S) {
    for (k in 1:D) {
      z_pred[s, k] = normal_rng(z[s, k], obs_sigma);
    }
  }
}
