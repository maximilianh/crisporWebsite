data {
  int<lower=1> N; // Number of data points 
  int<lower=18> L; // Length of the guide sequence
  int<lower=1> num_guides; // Number of guides present (equivalent to maximum value in guide_index)
  int<lower=1> matrix_cols; // L * mismatch_types

  int guide_index[N]; //Links observations to which guide they are from

  vector<lower=0, upper=1.0>[N] indel_rate; // Observed indel rate
  matrix[N, matrix_cols] mm_matrix; // Describes the occurrence of mismatches

  real beta_hyperparameter; # Hyperparameter that sets the regularization for the penalty parameters

  real<lower=0.0001, upper=10.0> sigma; #Standard deviation / noise of the indel rate

}

parameters {
  vector<lower=0>[matrix_cols] mm_penalty; 
  vector<lower=0, upper=1.0>[num_guides] guide_effect; 
}

transformed parameters {

  vector[N] total_penalty;
  vector<lower=0, upper=1.0>[N] pred_indel;
  vector[matrix_cols] mult_effect;

  mult_effect <- exp(-mm_penalty);
  total_penalty <- exp(-mm_matrix * mm_penalty);
  for (i in 1:N)  
    pred_indel[i] <- guide_effect[guide_index[i]] * total_penalty[i];

}

model {

  for (i in 1:matrix_cols)
    mm_penalty[i] ~ exponential(beta_hyperparameter);

  indel_rate ~ normal(pred_indel, sigma);

}

