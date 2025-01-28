// Started 27 January 2025
// By Ken

// Logistic model for damselfly speed as a function of wetmass

data{
  int <lower = 0> N;    // number of observations
  int <lower = 0> N_sp; // number of species
  
  real <lower = 0> wm[N];  // wetmass
  int <lower = 0> sp[N];  // wetmass
  real <lower = 0> spd[N]; // speed
}

parameters{
  real maxspd[N_sp];
  real mmax;
  real <lower = 0> smax;
  
  real k[N_sp];
  real mk;
  real <lower = 0> sk;
  
  real mp[N_sp];
  real mmp;
  real <lower = 0> smp;
  
  real <lower = 0> error;
}

transformed parameters{
  real spd_pred[N];
  for(i in 1:N){
    spd_pred[i] = maxspd[sp[i]] / (1 + exp(-k[sp[i]] * (wm[i] - mp[sp[i]])));
  }
}

model{
  maxspd ~ normal(mmax, smax);
  mmax ~ normal(30, 3 / 2.32);
  smax ~ normal(3 / 2.32, 3 / pow(2.32, 2));
  
  k ~ normal(mk, sk);
  mk ~ normal(0.2, 0.05 / 2.32);
  sk ~ normal(0.05 / 2.32, 0.05 / pow(2.32, 2));
  
  mp ~ normal(mmp, smp);
  mmp ~ normal(13, 3 / 2.32);
  smp ~ normal(4 / 2.32, 3 / pow(2.32, 2));
  
  error ~ normal(5, 3 / 2.32);
  
  for(i in 1:N){
    target += normal_lpdf(spd[i] | spd_pred[i], error);
  }
}
