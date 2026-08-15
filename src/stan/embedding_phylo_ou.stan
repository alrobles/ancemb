// embedding_phylo_ou.stan
// Bayesian phylogenetic model for protein embedding evolution
// using an Ornstein-Uhlenbeck process.
//
// Extends embedding_phylo_bm.stan by adding selection strength alpha
// and diffusion variance sigma2 per dimension. The process starts at a
// fixed root state z_anc, which is also the optimum (theta = z_anc).
//
// Fixed-root OU covariance between tips i, j for dimension k:
//   V[i,j,k] = (sigma2[k] / (2*alpha[k])) *
//              (exp(-alpha[k] * (t_i + t_j - 2*t_mrca(i,j)))
//               - exp(-alpha[k] * (t_i + t_j)))
//
// where t_i is the root-to-tip distance for species i.
//
// Non-centered parameterization uses the Cholesky of V_ou.

data {
  int<lower=1> S;
  int<lower=1> D;

  matrix[S, D] z_obs;
  matrix[S, S] C;                   // phylogenetic BM VCV (from ape::vcv)
  vector[S] tip_heights;            // root-to-tip distances

  real<lower=0> obs_sigma;

  vector[D] z_anc_prior_mean;
  vector<lower=0>[D] z_anc_prior_sd;
  real<lower=0> rate_scale;
  real<lower=0> alpha_scale;
}

parameters {
  vector[D] z_anc;
  vector<lower=0>[D] sigma2;
  vector<lower=0>[D] alpha;
  matrix[S, D] z_raw;
}

transformed parameters {
  matrix[S, D] z;

  for (k in 1:D) {
    matrix[S, S] V_ou;
    real halflife_inv = 2 * alpha[k] + 1e-10;
    real var_scale = sigma2[k] / halflife_inv;

    for (i in 1:S) {
      for (j in 1:S) {
        real shared = C[i, j];
        real unshared_i = tip_heights[i] - shared;
        real unshared_j = tip_heights[j] - shared;
        real total = tip_heights[i] + tip_heights[j];
        V_ou[i, j] = var_scale * (exp(-alpha[k] * (unshared_i + unshared_j))
                                  - exp(-alpha[k] * total));
      }
    }

    {
      matrix[S, S] V_ou_jitter = V_ou + diag_matrix(rep_vector(1e-8, S));
      matrix[S, S] L_ou = cholesky_decompose(V_ou_jitter);
      z[, k] = z_anc[k] + L_ou * z_raw[, k];
    }
  }
}

model {
  z_anc ~ normal(z_anc_prior_mean, z_anc_prior_sd);
  sigma2 ~ normal(0, rate_scale);
  alpha ~ normal(0, alpha_scale);

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
