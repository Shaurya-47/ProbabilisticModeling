library(rstan)
rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())
library(xts)
library(coda)
library(rjags)
library(outbreaks)
library("lubridate")
library(ggplot2)

data1 <- ebola_kikwit_1995
data1$Day <- yday(data1$date)

#Daily case count plot
ggplot(data1, aes(x=Day, y=onset)) + 
  geom_bar(stat = "identity", width=0.5) +
  theme_classic()+ scale_x_continuous(name="Time of symptom onset (days)",expand = c(0, 0), limits = c(0,NA), breaks = seq(0,200,20)) + 
  scale_y_continuous(name= "Cases",expand = c(0, 0), limits = c(0, NA), breaks = seq(0,20,5))
#Daily death count plot
ggplot(data1, aes(x=Day, y=death)) + 
  geom_bar(stat = "identity", width=0.5) +
  theme_classic()+ scale_x_continuous(name="Time of death (days)",expand = c(0, 0), limits = c(0,NA), breaks = seq(0,200,20)) + 
  scale_y_continuous(name= "Cases",expand = c(0, 0), limits = c(0,15),breaks = seq(0,15,5))

C<- data1$onset
D<- data1$death
n_days <- length(C)
ncases <- sum(C)
s0=5364500-1  #initial conditions
e0=1          #initial conditions
a=0           #initial conditions
beta=0.2      #initial conditions
q=0.2         #initial conditions
rho=0.2       #initial conditions
gamma=0.143   #initial conditions
h=1           #initial conditions
r0=0          #initial conditions
tcontrol=130  # time point when control interventions was introduced
N = 5364500   # population size

ST= N-ncases
m= ST-s0
abs(m) 
S <- vector(mode="numeric", length = 192) #Susceptible compartment
S[1]=s0 #  initial conditions for S
E <- vector(mode="numeric", length = 192) #Exposed compartment
E[1]= e0 #  initial conditions for E
I <- vector(mode="numeric", length = 192) #Infectious compartment
I[1]= a #  initial conditions for I
R <- vector(mode="numeric", length = 192) #Recovered(Dead or recovered) compartment
R[1]= r0 #  initial conditions for R
B0 <- vector(mode="numeric", length = 192) #Number of susceptible who become infectious
#y0 = c(S = s0,E=e0, I = a, R = r0)
B0[1] <-abs(m)
Beta <- vector(mode="numeric", length = 192)
B <- vector(mode="numeric", length = 192)
P <-vector(mode="numeric", length = 192)
pc<-vector(mode="numeric", length = 192)
pr <-vector(mode="numeric", length = 192)

B[1] = B0[1]
S[1]= S[1] - B0[1]
E[1]= E[1] + B0[1] - C[1]
I[1]= I[1] + C[1] - D[1]
R[1]= N - S[1] - E[1] - I[1]
# Simulating B using initial condition and empirical data of C and D 
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
  
  B[t] = B[t] = rbinom(n=1,size=S[t], prob=P[t])
}

#input for the Stan code
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
  real<lower=0,upper=1> beta; //parameter constraints mentioned in the paper are used as bounds
  real<lower=1/10.7,upper=1/3.5> gamma;
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
  beta ~ gamma(2,10); // prior
  gamma ~  gamma(2,14); // prior
  rho ~  gamma(2,10); // prior
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
      B[t] ~ binomial(S[t], P[t]); // # of cases infected but not yet infectious.
      C[t] ~ binomial(E[t], pc); // # of cases by date of symptom onset
      D[t] ~ binomial(I[t], pr); // # of cases removed from infectious class
  }
}
generated quantities {
  real Reproduction = beta/gamma; // reproduction number
  real EffReproduction[n_days];   // effective reproduction rate
  for(t in 130:n_days){
  EffReproduction[t] = beta*exp(-q*(t-130))/gamma; 
  }
}",
"stan_model1.stan")
stanc("stan_model1.stan")
stan_model1 <- "stan_model1.stan"
parameter= c("beta", "gamma", "rho", "q", "Reproduction", "EffReproduction")

fit_seir <- stan(file=stan_model1,
                 data = DataSEIR,
                 pars = parameter,
                 #init = ini,
                 iter = 20000,
                 chains = 2,
                 warmup = 10000)
summary(fit_seir, pars=c("beta", "gamma", "rho", "q", "Reproduction"))
stan_dens(fit_seir, pars=c("beta", "gamma", "rho", "q", "Reproduction"))


a <- summary(fit_seir, pars=c("EffReproduction"))
b = a[["summary"]][,c(1,3)] # extracting the mean and sd
plot(b[131:192,1],xlab="Time",ylab="Values",main="Effective R0")
