require(rstan)
require(MASS)
library(parallel)
require(beeswarm)
setwd("/home/yi/Dropbox/EvolutionOfVariance/CodeForSubmission/")

#Simulation
##The parameter names follow the notation in the main text
n.ind = 100 ##number of individuals
n.reps<-20 ##repeated measures per individual

Bd<-0   ## developmental plasticity
Br<-1   ## reversible plasticity

Va<-3   ## among individual variance in the means
VBd<-0  ## variance in developmental plasticity
VBr<-1  ## variance in reversible plasticity, also variance in within individual standard deviation and within individual variance. 

beta<-0## directional selection
delta<-  0.15 ## non-linear selection
eta <-  -0.35   ## selection on within-individual variance

Vxd<-0  ## variance in the developmental environment
Vr<-7   ## residual variance
Vxr<-Vr/(VBr + Br^2) ## Variance in the transient environment/current environment.

## P matrix
m<-matrix(0,3,3)
m[1,1]<-Va
m[2,2]<-VBd
m[3,3]<-VBr

##Create data set with individual specific values  
d<-as.data.frame(mvrnorm(n.ind, c(0, Bd, Br), m))
colnames(d)<-c("a", "Bd", "Br")
d$ID<-1:n.ind
d$xd<-rnorm(nrow(d),0,sqrt(Vxd))
d$i<-d$a + d$Bd*d$xd
d$v <-  beta*d$a +  delta*d$a^2 + eta*d$Br^2*Vxr #+ sf 
d$w<-rpois(nrow(d),exp(d$v))

##Create the data set with repeated measures  
d2<-d[rep(seq_len(nrow(d)), n.reps), ]
d2$xr<-rnorm(nrow(d2),0,sqrt(Vxr))
d2$z<-d2$i + d2$Br*d2$xr

#Using Stan to fit the model
##Write Stan model within R  
write(
  "data {
    // Define variables in data
     int<lower=0> N; //number of observation
     int<lower=0> nind; //number of individuals
     
     // Cluster IDs
     int<lower=0> ID[N]; // Individual IDs
     
     // Continuous outcomes
      real t[N]; //Phenotypic observations
      int<lower=0> w[nind]; //Fitness measure  
    }
    
    parameters {
      // Define parameters to estimate
      //For the phenotype
      //Fixed effects
     real c;  // Population intercept
     real m_lsigma_e; //Average log residual standard deviation
    // Random effects
       matrix[2,nind] zI; // intercepts and log residual standard deviation for each individual
    
     //For fitness
     real alpha; //Mean log fitness
     real beta; // selection on individual means
     real delta; // selection on individual squared deviations from the population mean
     real eta; // selection on an individual's squared deviations from its own mean
      
      // Random effects
       vector<lower=0>[2]      sigma_I; // intercepts and log residual standard deviation for each individual
       cholesky_factor_corr[2] L;  // factor to estimate covariance between intercepts and log residual standard deviation for each individual
    }
    
   transformed parameters{
    real<lower=0> sigma_e[nind];
    row_vector[nind] sigma2_e;
    real mu_t[N];
    real mu_v[nind];
      
    matrix[2,nind] I; //  Unsclaed blups intercept and slope for w island year
    I  = diag_pre_multiply(sigma_I, L) * zI; // get the unscaled value
    
   for (i in 1:nind){
     sigma_e[i] = exp(m_lsigma_e + I[2, i]);
     sigma2_e[i] = sigma_e[i]^2;
      mu_t[i] = c + I[1,i];
      mu_v[i] = alpha + beta*I[1,i] + delta*I[1,i]^2 + eta*sigma2_e[i];
                    }
                      }

   model{
      // Random effects distribution
      to_vector(zI) ~ normal(0, 1);
      alpha ~ normal(0, 1);
      beta ~ normal(0, 1);
      delta ~ normal(0, 1);
      eta ~ normal(0, 1);
      
      to_vector(sigma_I) ~ normal(0, 1);
      m_lsigma_e ~ normal(0, 1);
      L ~ lkj_corr_cholesky(2);
  
      // Likelihood part of Bayesian inference
      for (i in 1:N) {
        t[i] ~ normal(mu_t[ID[i]], sigma_e[ID[i]]);    
                     }
    
        w ~ poisson(exp(mu_v));    
                      
                       }"
  , "ME2.stan")


##Store data sets in a list
stan_data_list <- list(t = d2$z,  nind = n.ind, 
                       w= d$w, 
                       N=nrow(d2), ID=d2$ID, v=d$v)

## MCMC settings
ni <- 15000
nt <- 10
nb <- 5000
nc <- 3

##Parameters monitored
params <- c("beta", "delta", "eta") ##These are the different selection gradients

##Fitting the model
md <- stan("ME2.stan", data = stan_data_list,  pars = params,
           chains = nc, iter = ni, warmup = nb, thin = nt, cores=3) ##"ME2.stan" is name and path for the stan model

summary(md)$summary
