library(AGHmatrix)
library(hglm)
library(rstan)
setwd("/home/yi/Dropbox/Evolutionofvariance")

Br=1
Bd=1
xd_bar=0
xr_bar=0

Va=2

VBd=0
Vxd=1

VBr=0
Vxr=1

res1<-simPed(n=200, gen=4, nrep=4, Va=Va, VBd=VBd, Vxd=Vxd, VBr=VBr, Vxr=Vxr, Br=Br, Bd=Bd, xd_bar=xd_bar, xr_bar=xr_bar)
x<-res1$d2

ped<-x[-which(duplicated(x$ID)),c(3,2,1)]
#ped<-x[,c(3,2,1)]
colnames(ped)<-c("animal","dam", "sire")
#ped$sire[ped$sire==0]<-NA
#ped$dam[ped$dam==0]<-NA  
#colnames(ped)<-c("ID","dam.id","sire.id")
library(tidyr)

A  <- Amatrix(ped)
Z0  <- diag(nrow(ped))
L <- t(chol(A))
Z  <- Z0 %*% L

stan_data = list(Y = x$z,
            A = A,
            Z = Z0,
            K = max(x$ID),
            N = nrow(x),
            ID=x$ID)

inits <- function() list(sigma_G = runif(1, 0.01, 10),
                         sigma_R = runif(1, 0.01, 10),
                         sigma_I = runif(1, 0.01, 10))

ni <- 1000
nt <- 1
nb <- 500
nc <- 1
params<-c("sigma_U", "sigma_E")

## Call Stan from R
md <- stan("ME1.stan", data = stan_data, init = inits, pars = params,
           chains = nc, iter = ni, warmup = nb, thin = nt, cores=1)


summary(md)
summary(md)$summary



write(
  "data {
  int<lower=1>    K; // number of all animals
  int<lower=1>    N; // number of observations
  matrix[N,K]     Z; // Random effects design matrix
  vector[N]       Y; // response variable
  matrix[K,K]     A; // relationship matrix
}
transformed data{
  matrix[K,K] LA;
  LA = cholesky_decompose(A);
}
parameters {
  vector[K]  a_decompose; // breeding values
  real<lower=0> sigma_G; // genetic standard deviation
  real<lower=0> sigma_R; // residual standard deviation
}
model {
  vector[K] muI;
  vector[K] a;

  a_decompose ~ normal(0, 1);
  
  a = sigma_G * (LA * a_decompose);
  muI =  Z * a + Z*;
  
  Y~normal(muI, sigma_R);
  
  sigma_G ~ inv_gamma_lpdf( sigma_G  | 1, 1);
  sigma_R ~ inv_gamma_lpdf( sigma_G  | 1, 1);
}
generated quantities{
  real sigma_U;
  real sigma_E;

  sigma_U = sigma_G * sigma_G; // genetic variance
  sigma_E = sigma_R * sigma_R; // residual variance
}"
, "ME1.stan")


Br=1
Bd=1
xd_bar=0
xr_bar=0

Va=1

VBd=0
Vxd=1

VBr=0
Vxr=1

res1<-simPed(n=50, gen=4, nrep=5, Va=Va, VBd=VBd, Vxd=Vxd, VBr=VBr, Vxr=Vxr, Br=Br, Bd=Bd, xd_bar=xd_bar, xr_bar=xr_bar)
x<-res1$d2
x$ob<-1:nrow(x)

Z<-table(x$ob, x$ID)

ped<-x[-which(duplicated(x$ID)),c(3,2,1)]
colnames(ped)<-c("animal","dam", "sire")

A  <- Amatrix(ped)

stan_data = list(Y = x$z,
                 A = A,
                 Z = Z,
                 K = max(x$ID),
                 N = nrow(x),
                 ID=x$ID)

inits <- function() list(sigma_G = runif(1, 0.01, 2),
                         sigma_R = runif(1, 0.01, 2),
                         sigma_I = runif(1, 0.01, 2))

ni <- 1000
nt <- 2
nb <- 500
nc <- 2
params<-c("sigma2_I", "sigma2_E", "sigma2_G")

## Call Stan from R
md2 <- stan("ME2.stan", data = stan_data, init = inits, pars = params,
           chains = nc, iter = ni, warmup = nb, thin = nt, cores=3)


summary(md2)$summary


write(
  "data {
  int<lower=1>    K; // number of all animals
  int<lower=1>    N; // number of observations
  matrix[N,K]     Z; // Random effects design matrix
  vector[N]       Y; // response variable
  matrix[K,K]     A; // relationship matrix
  int<lower=1>    ID[N]; 
}
transformed data{
  matrix[K,K] LA;
  LA = cholesky_decompose(A);
}
parameters {
  vector[K]  a_decompose; // breeding values
  vector[K]  I_z; // breeding values
  
  real<lower=0> sigma_G; // genetic standard deviation
  real<lower=0> sigma_R; // residual standard deviation
  real<lower=0> sigma_I; // PE standard deviation
}
model {
  vector[N] muI;
  vector[K] a;
  vector[K] I;
  
  a_decompose ~ normal(0, 1);
  I_z ~ normal(0, 1);
  
  a = sigma_G * (LA * a_decompose);
  I = sigma_I * I_z;
  muI =  Z*a + Z*I;
  
  
target += normal_lpdf(Y  | muI, sigma_R);

target +=inv_gamma_lpdf( sigma_I  | 1, 1);
target +=inv_gamma_lpdf( sigma_R  | 1, 1);
target += inv_gamma_lpdf( sigma_G  | 1, 1);
}

generated quantities{
  real sigma2_G;
  real sigma2_E;
  real sigma2_I;
  sigma2_G = sigma_G * sigma_G; // genetic variance
  sigma2_I = sigma_I * sigma_I; // genetic variance
  sigma2_E = sigma_R * sigma_R; // residual variance
}"
, "ME2.stan")


Br=1
Bd=0
xd_bar=0
xr_bar=0

Va=2

VBd=0
Vxd=0

VBr=0
Vxr=0.5

res1<-simPed(n=50, gen=5, nrep=6, Va=Va, VBd=VBd, Vxd=Vxd, VBr=VBr, Vxr=Vxr, Br=Br, Bd=Bd, xd_bar=xd_bar, xr_bar=xr_bar)
x<-res1$d2
x$ob<-1:nrow(x)

ped<-x[-which(duplicated(x$ID)),c(3,2,1)]
colnames(ped)<-c("animal","dam", "sire")

A  <- Amatrix(ped)
Z<-table(x$ob, x$ID)
Z2 <- diag(nrow(ped))

stan_data = list(Y = x$z,
                 Y2 = x$z^2,
                 A = A,
                 Z = Z,
                 Z2 = Z2,
                 K = max(x$ID),
                 N = nrow(x),
                 ID=x$ID)

inits <- function() list(sigma_G = runif(1, 0.01, 3),
                         sigma_R = runif(1, 0.01, 3),
                         sigma_I = runif(1, 0.01, 3))

ni <- 2000
nt <- 2
nb <- 1000
nc <- 1
params<-c("sigma2_G", "sigma2_I", "sigma2_E")

## Call Stan from R
md3 <- stan("/home/yi/Dropbox/Evolutionofvariance/ME3.stan", data = stan_data, init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt, cores=1)


summary(md3)$summary


xwrite(
  "data {
  int<lower=1>    K; // number of individuals
  int<lower=1>    N; // number of observations
  matrix[N,K]     Z; // Random effects design matrix
  matrix[K,K]     Z2; // Random effects design matrix2
  vector[N]       Y; // response variable
  vector[N]       Y2; // response variable
  matrix[K,K]     A; // relationship matrix
}

transformed data{
  matrix[K,K] LA;
  LA = cholesky_decompose(A);
}

parameters {
  vector[K]  a_decompose; // standardized breeding values
  vector[K] I; // individual deviations 
  vector[K] I2b; // individual deviations 
  
  real b1; //mean of observed trait
  real b2; //mean of observed trait

  real<lower=0> sigma_I; // Individual standard deviation
  real<lower=0> sigma_R; // residual standard deviation
  real<lower=0> sigma_R2; // residual standard deviation
 
  real<lower=0> sigma_G; // genetic standard deviation for squared values
  real<lower=0> sigma_I2; // environment standard deviation for sq values
}

transformed parameters{
  vector[K] I2; // squared of individual deviations
  for(j in 1:K){ // Calculate squared values of the mean deviations
  I2[j] =  I[j];  
   }
}

model {
vector[N] mu1; // expected observation values
vector[N] mu2; // expected squared values
vector[K] muI2; // expected mean squared values
vector[K] a;  // breeding values for squared individual values

a = sigma_G * (LA * a_decompose);
target += normal_lpdf(a_decompose  | 0, 1); 
target += normal_lpdf(I  | 0, sigma_I);   

mu1 = b1 + Z*I; // Effect of the individual average on phenotype expression
muI2 =  b2 + Z2*a;  // Effect of a on the average squared deviation
mu2 = Z*muI2;
  
target += normal_lpdf(Y  | mu1, sigma_R);   //model for phenotype
target += normal_lpdf(I2  | muI2, sigma_I2); //model for squared
target += normal_lpdf(Y2  | mu2, sigma_R2);  //model for squared
 
//Priors
b1~normal(0,1);
b2~normal(0,1);

target +=inv_gamma_lpdf( sigma_I  | 1, 1);
target +=inv_gamma_lpdf( sigma_R  | 1, 1);
target +=inv_gamma_lpdf( sigma_R2  | 1, 1);

target += inv_gamma_lpdf( sigma_G  | 1, 1);
target +=inv_gamma_lpdf( sigma_I2  | 1, 1);
 
}

generated quantities{
  real sigma2_G;
  real sigma2_E; 
  real sigma2_I;
  
  sigma2_G = sigma_G*sigma_G;
  sigma2_I = sigma_I * sigma_I; // individual variance
  sigma2_E = sigma_R * sigma_R; // residual variance
}"
, "/home/yi/Dropbox/Evolutionofvariance/ME3.stan")



write(
  "data {
  int<lower=1>    K; // number of individuals
  int<lower=1>    N; // number of observations
  matrix[N,K]     Z; // Random effects design matrix
  matrix[K,K]     Z2; // Random effects design matrix
  vector[N]       Y; // response variable
  matrix[K,K]     A; // relationship matrix
}

transformed data{
  matrix[K,K] LA;
  LA = cholesky_decompose(A);
}

parameters {
  vector[K]  a_decompose; // standardized breeding values
  vector[K] I; // individual deviations 
  
  real b1; //mean of observed trait
  real b2; //mean of observed trait

  real<lower=0> sigma_I; // Individual standard deviation
  real<lower=0> sigma_R; // residual standard deviation
  
  real<lower=0> sigma_G; // genetic standard deviation for squared values
  real<lower=0> sigma_I2; // environment standard deviation for sq values
}

transformed parameters{
  vector[K] I2; // squared of individual deviations
  for(j in 1:K){ // Calculate squared values of the mean deviations
  I2[j] =  I[j]^2;  
   }
}

model {
vector[N] mu1; // expected observation values
vector[K] muI2; // expected mean squared values
vector[K] a;  // breeding values for squared individual values

a = sigma_G * (LA * a_decompose);
target += normal_lpdf(a_decompose  | 0, 1); 
target += normal_lpdf(I  | 0, sigma_I);   
   
mu1 = b1 + Z*I; // Effect of the individual average on phenotype expression
muI2 =  b2 + Z2*a;  // Effect of a on the average squared deviation

target += normal_lpdf(Y  | mu1, sigma_R);   //model for phenotype
target += normal_lpdf(I2  | muI2, sigma_I2); //model for squared

//Priors
b1~normal(0,1);
b2~normal(0,1);

target +=inv_gamma_lpdf( sigma_I  | 1, 1);
target +=inv_gamma_lpdf( sigma_R  | 1, 1);

target += inv_gamma_lpdf( sigma_G  | 1, 1);
target +=inv_gamma_lpdf( sigma_I2  | 1, 1);
 
}

generated quantities{
  real sigma2_G;
  real sigma2_E; 
  real sigma2_I;
  
  sigma2_I = sigma_I * sigma_I; // individual variance
  sigma2_E = sigma_R * sigma_R; // residual variance
}"
, "/home/yi/Dropbox/Evolutionofvariance/ME3.stan")





