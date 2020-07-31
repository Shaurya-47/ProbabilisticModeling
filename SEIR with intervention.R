# SEIR with intervention

library(rstan)
rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')
library(outbreaks)


data1 <- ebola_kikwit_1995
C<- data1$onset
D<- data1$death
n_days <- length(C)
ncases <- sum(C)
s0=5364500-1
e0=1 
a=0
beta=0.2
q=0.2
rho=0.2
gamma=0.143 
h=1
r0=0
tcontrol=130
N = 5364500
t <- seq(0, n_days, by = 1)
t0 = 0 
t <- t[-1]
ST= N-ncases
m= ST-s0
abs(m) #  final size of the epidemic

S <- vector(mode="numeric", length = 192)
S[1]=s0 #  initial conditions for S
E <- vector(mode="numeric", length = 192)
E[1]= e0 #  initial conditions for E
I <- vector(mode="numeric", length = 192)
I[1]= a #  initial conditions for I
R <- vector(mode="numeric", length = 192)
R[1]= r0 #  initial conditions for R

B0 <- vector(mode="numeric", length = 192)

B0[1] <- abs(m)
Beta <- vector(mode="numeric", length = 192)
B <- vector(mode="numeric", length = 192)
P <-vector(mode="numeric", length = 192)
pc<-vector(mode="numeric", length = 192)
pr <-vector(mode="numeric", length = 192)
B[1] = 0 
S[1]= S[1] - B0[1]
E[1]= E[1] + B0[1] - C[1]
I[1]= I[1] + C[1] - D[1]
R[1]= N - S[1] - E[1] - I[1]
for(t in 2:n_days){
  if(t >= 130){
    Beta[t]= beta*exp(-q*(t-130))
  }
  else{
    Beta[t]= beta
  }
  
  S[t]= S[t-1] - B[t-1]
  E[t]= E[t-1] + B[t-1] - C[t-1]
  I[t]= I[t-1] + C[t-1] - D[t-1]
  R[t] = N - S[t] - E[t] - I[t]
  
  P[t]= 1 - exp(-(Beta[t]*I[t]/N))
  pc = 1 - exp(-rho)
  pr = 1- exp(-gamma)
  
  B[t] = rbinom(n=1,size=S[t], prob=P[t])
}

DataSEIR <- list(n_days = n_days, S0 = S , E0 = E, I0 = I, R0 = R,
                 N = N, tcon = tcontrol, C = C,
                 D = D, B = B)

write("
data {
  int<lower=1> n_days;// total # of time point 
  int S0[n_days]; //of susceptible individuals at the start of the epidemic 
  int E0[n_days]; //of exposed individuals at the start of the epidemic
  int I0[n_days]; //of infectious individuals at the start of the epidemic
  int R0[n_days]; //of removed individuals at the start of the epidemic
  int N;          //population size
  int C[n_days];  // empirical data
  int D[n_days];  // empirical data
  int tcon;       // time point where control measures are introduced 
  int B[n_days]; // initial conditions for unobserved entity. 
}
parameters {
  real<lower=0,upper=1> beta;
  real<lower=1/10.5,upper=1/3.5> gamma;
  real<lower=1/21.0,upper=1> rho;
  real<lower=0> q;
}
transformed parameters{
  real theta[4];
  theta[1]= beta;
  theta[2]= gamma;
  theta[3]= rho;
  theta[4]= q;
}
model{
  int S[n_days]; // # of susceptible individuals
  int E[n_days]; // # of exposed individuals
  int I[n_days]; // # of infectious individuals
  int R[n_days]; // # of removed individuals
  real P[n_days]; //  probability of leaving S compartment
  
  real Beta[n_days];
  real pc; //probability of leaving E compartment
  real pr; //probability of leaving I compartment
  
  
  S[1]= S0[1] - B[1];
  E[1]= E0[1] + B[1] - C[1];
  I[1]= I0[1] + C[1] - D[1];
  R[1]= N - S[1] - E[1] - I[1];
  
  beta ~ gamma(20,100); // prior
  gamma ~  gamma(20,140); // prior
  rho ~  gamma(20,100); // prior
  q ~ gamma(2,10); // prior
  
  for(t in 2:n_days){
      if(t >= tcon){
        Beta[t]= beta*exp(-q*(t-tcon));
      }
      else{
        Beta[t]= beta;
      }
      
      S[t]= S[t-1] - B[t-1];
      E[t]= E[t-1] + B[t-1] - C[t-1];
      I[t]= I[t-1] + C[t-1] - D[t-1];
      R[t] = N - S[t] - E[t] - I[t];
      
      P[t]= 1 - exp(-(Beta[t]*I[t]/N));
      pc = 1 - exp(-rho);
      pr = 1- exp(-gamma);
      
      
      B[t] ~ binomial(S[t], P[t]); //  
      C[t] ~ binomial(E[t], pc); // # of cases by date of symptom onset
      D[t] ~ binomial(I[t], pr); // # of cases removed from infectious class
      
      
  }
}

generated quantities {
  real Reproduction = beta/gamma; // reproduction number
  real EffReproduction[n_days];   // effective reproduction rate
  for(t in 131:n_days){
  EffReproduction[t] = beta*exp(-q*(t-tcon))/gamma; 
  }
}
",
"stan_model1.stan")
stanc("stan_model1.stan")
stan_model1 <- "stan_model1.stan"

parameter= c("beta","gamma","rho","q", "Reproduction", "EffReproduction")
#ini= function(){list(beta = 0.2, gamma = 0.143, rho=0.2, q = 0.2)}
fit_seir <- stan(file=stan_model1,
                 data = DataSEIR,
                 pars = parameter,
                 #init = ini,
                 iter = 20000,
                 chains = 2,
                 warmup = 10000)


# Summary object for only the three parameters and Basic Reproduction number
summary(fit_seir, pars=c("beta", "gamma", "rho", "q", "Reproduction"), probs = c(0.1,0.9))
stan_dens(fit_seir, pars=c("beta", "gamma", "rho", "q", "Reproduction"))

# Summary statistics for the Effective Reproduction number - varies over time
a <- summary(fit_seir, pars=c("EffReproduction"))
b = a[["summary"]][,c(1,3)] # extracting the mean and sd
plot(b[131:192,1], main="Effective Reproduction Number", xlab="time", ylab="value")
