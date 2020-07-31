# SIR without intervention

#importing libraries
library(rstan)
rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')
library(outbreaks)

#data and initialisations
data1 <- ebola_kikwit_1995
C <- data1$onset 
D <- data1$death 
n_days <- length(C)
ncases <- sum(C)
s0=  5364500-1 
i0=1
r0=0
h=1
N = 5364500
t <- seq(0, n_days, by = 1)
t0 = 0 
t <- t[-1]
ST= N-ncases
m= ST-s0
abs(m) #  final size of the epidemic

S <- vector(mode="numeric", length = 192) 
S[1]=s0 #  initial conditions for S
I <- vector(mode="numeric", length = 192) 
I[1]= i0 #  initial conditions for I
R <- vector(mode="numeric", length = 192) 
R[1]= r0 #  initial conditions for R

# data for Stan
DataSEIR <- list(n_days = n_days, S0 = S, I0 = I, R0 = R,
                 N = N, C = C, D = D)

# Stan code
write("
data {
  int<lower=1> n_days;// total number of of time points
  int S0[n_days]; //of susceptible individuals at the start of the epidemic 
  int I0[n_days]; //of infectious individuals at the start of the epidemic
  int R0[n_days]; //of removed individuals at the start of the epidemic
  int N;          //population size
  int C[n_days];  // empirical data
  int D[n_days];  // empirical data
}
parameters {
  real<lower=0,upper=1> beta;
  real<lower=1/10.7,upper=1/3.5> gamma;
}
transformed parameters{
  real theta[2];
  theta[1]= beta;
  theta[2]= gamma;
}
model{
  int S[n_days]; // # of susceptible individuals
  int I[n_days]; // # of infectious individuals
  int R[n_days]; // # of removed individuals

  S[1]= S0[1] - C[1];
  I[1]= I0[1] + C[1] - D[1];
  R[1]= N - S[1] - I[1];
  
   beta ~ gamma(20,100); // prior
   gamma ~  gamma(20,140); // prior
  
  for(t in 2:n_days){
      S[t]= S[t-1] - C[t-1];
      I[t]= I[t-1] + C[t-1] - D[t-1];
      R[t] = N - S[t] - I[t];
      
      C[t] ~ binomial(S[t],1 - exp(-beta*I[t]/N)); 
      D[t] ~ binomial(I[t],1- exp(-gamma)); 
    }

}
generated quantities {
  real Reproduction = beta/gamma;
}

", "stan_model1.stan")
stanc("stan_model1.stan")
stan_model1 <- "stan_model1.stan"

# fitting the model
parameter= c("beta","gamma","Reproduction")
#ini= function(){list(beta = 0.2, gamma = 0.143)}
fit_sir <- stan(file=stan_model1,
                 data = DataSEIR,
                 pars = parameter,
                 #init = ini,
                 iter = 20000,
                 chains = 2, warmup=10000)

#summary
summary(fit_sir, pars=c("beta", "gamma", "Reproduction"), probs = c(0.1,0.9))
stan_dens(fit_sir, pars=c("beta", "gamma", "Reproduction"))
