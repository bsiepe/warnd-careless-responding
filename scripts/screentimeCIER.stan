data{
  int<lower = 1> N; // number of persons
  int<lower = 1> Nobs; // number of observations
  int<lower = 1> person[Nobs]; // person id
  int<lower = 1> obs_id[Nobs]; // number of prompts responded to since the start of the EMA
  int<lower = 1> prompt_num[Nobs]; // number of prompts received since the start of the day
  int<lower = 1> day_num[Nobs]; // number of days passed since start of the EMA
  real log_total_time[Nobs]; // log total time on the target segment
}
parameters{
  row_vector[3] PersPar[N]; // person parameters: 1: attentiveness, 2: log initial time expenditure, 3: log habitual decay
  real gamma0; // initial attentiveness difficulty
  real gammaD; // change of attentiveness difficulty across days
  real gammaB; // change of attentiveness difficulty within days (across prompts)
  corr_matrix[3] correlP; // correlations of person parameters
  real muC; // average careless log time
  real<lower=0> offsetA; // offset, capturing difference between average careless log time and lower bound of attentive times
  real muTau1; // mean of log initial time expenditure
  real muTau0; // mean of log habitual decay
  real<lower=0> sigmaC; // standard deviation of careless times
  real<lower=0> sigmaA; // residual standard deviation of attentive times
  vector<lower=0>[3] sigmaP; // sd person parameters
}
transformed parameters{
  real<lower=0,upper=1> pDelta[Nobs]; // attentiveness probabilities  
  real eta[Nobs]; // argument for attentiveness probabilities logit  
  cov_matrix[3] SigmaP; // person parameter variance covariance matrix
  SigmaP = quad_form_diag(correlP, sigmaP);
  for(n in 1:Nobs){
    eta[n]=PersPar[person[n],1] -(gamma0 + gammaD*(day_num[n]-1) + gammaB*(prompt_num[n]-1));
    pDelta[n] = 1/(1+exp(-eta[n])); 
  }
}
model{
  // priors
  PersPar~ multi_normal(rep_vector(0,3),SigmaP); 
  correlP ~ lkj_corr(1); 
  sigmaP ~  cauchy(0,2); 
  muTau0  ~  normal(0,5); 
  muTau1  ~  normal(0,5);   
  
  gamma0 ~ normal(0, 5);
  gammaD ~ normal(0, 5);
  gammaB ~ normal(0, 5);
  
  sigmaC ~  cauchy(0,2); 
  sigmaA ~  cauchy(0,2); 

  muC ~  normal(0,5); 
  offsetA ~  normal(0,5); 

  // model
  for(n in 1:Nobs) {
    target += log_mix(pDelta[n],
    normal_lpdf(log_total_time[n]|muC + offsetA + exp(PersPar[person[n],2]+muTau0)*(exp(PersPar[person[n],3]+muTau1))^(obs_id[n]-1), sigmaA), // attentive observation
    normal_lpdf(log_total_time[n]|muC, sigmaC));  // careless observation
  }
}
