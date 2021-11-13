data {
  int<lower=1>    K; // number of individuals
  int<lower=1>    N; // number of observations
  matrix[N,K]     Z; // Random effects design matrix
  matrix[K,K]     Z2; // Random effects design matrix2
  vector[N]       Y; // response variable
  matrix[K,K]     A; // relationship matrix
}

transformed data{
  matrix[K,K] LA;
  LA = cholesky_decompose(A);
}

parameters {
  vector[K]  a_decompose; // standardized breeding values
  vector[K]  I_z; // standardized breeding values
  

  real<lower=0> sigma_I; // Individual standard deviation
  real<lower=0> sigma_R; // residual standard deviation
  real<lower=0> sigma_G; // genetic standard deviation for squared values

  real b1; //mean of observed trait
}

transformed parameters{
vector[K] I;  // breeding values for squared individual values
vector[K] a;  // breeding values for squared individual values
a = sigma_G * (LA * a_decompose);
I = sigma_I * I_z;
}

model {
vector[N] mu1; // expected observation values

target += normal_lpdf(a_decompose  | 0, 1); 
target += normal_lpdf(I_z  | 0, 1); 

mu1 = b1 + Z*a + Z*I; // Effect of the individual average on phenotype expression

target += normal_lpdf(Y  | mu1, sigma_R);   //model for phenotype

//Priors
b1~normal(0,1);

target +=inv_gamma_lpdf(sigma_I  | 1, 1);
target +=inv_gamma_lpdf(sigma_R  | 1, 1);
target += inv_gamma_lpdf(sigma_G  | 1, 1);
}

generated quantities{
  real sigma2_G;
 real sigma2_I;
  real sigma2_R;
  
  sigma2_G = sigma_G*sigma_G;
  sigma2_I = sigma_I * sigma_I; // individual variance
  sigma2_R = sigma_R * sigma_R; // residual variance
}
