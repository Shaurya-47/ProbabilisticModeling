# Trying to construct the SIR model and simulate data from it (with intervention)

# Importing a package for differential equations - can be extended to discrete time as well
library(odin)
library(dde)
library(scales)

# i indicates which simulation
# Generating the Discrete Stochastic SEIR model without intervention dynamics
seir_mult_generator <- odin::odin({
  ## Core equations for transitions between compartments:
  update(S[]) <- S[i] - B[i]
  update(E[]) <- E[i] + B[i] - C[i]
  update(I[]) <- I[i] + C[i] - D[i]
  update(R[]) <- R[i] + D[i]
  
  ## Individual probabilities of transition: (intervention part added here)
  p_SE[] <- if (step>130) 1 - exp(-beta*exp(-q*(step-130))*I[i]/N[i]) else 1 - exp(-beta*I[i]/N[i]) # S to E
  p_EI <- 1 - exp(-rho)       # E to I
  p_IR <- 1 - exp(-gamma)     # I to R
  
  ## Draws from binomial distributions for numbers changing between compartments:
  B[] <- rbinom(S[i], p_SE[i])
  C[] <- rbinom(E[i], p_EI)
  D[] <- rbinom(I[i], p_IR)
  
  ## Total population size
  N[] <- S[i] + E[i] + I[i] + R[i]
  
  ## Initial states:
  initial(S[]) <- S_ini
  initial(E[]) <- E_ini
  initial(I[]) <- I_ini
  initial(R[]) <- R_ini
  
  ## User defined parameters - default in parentheses:
  S_ini <- user(5364499)
  E_ini <- user(1)
  I_ini <- user(0)
  R_ini <- user(0)
  beta <- user(0.17261446)
  rho <- user(0.04762486)
  gamma <- user(0.09374695)
  q <- user(0.18766678)
  
  ## Number of replicates
  nsim <- user(500)
  dim(N) <- nsim
  dim(S) <- nsim
  dim(E) <- nsim
  dim(I) <- nsim
  dim(R) <- nsim
  dim(p_SE) <- nsim
  dim(B) <- nsim
  dim(C) <- nsim
  dim(D) <- nsim
  
}, verbose = FALSE)


# Executing the simulation of a scenario
x <- seir_mult_generator()


# Generating a graph


set.seed(50)
seir_col <- c("gold", "darkorange", "maroon", "sky blue")
x_res_int <- x$run(0:500)
par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
matplot(x_res_int[, 1], x_res_int[, -1][,1:1500], xlab = "Time", ylab = "Number of individuals",
        type = "l", col = alpha(rep(seir_col, each = 500), 0.5), lty = 1, 
        mgp=c(3.2,0.4,0), ylim = c(0,200), xlim= c(0,400))
legend("topright", lwd = 1, col = seir_col, legend = c("S", "E", "I", "R"), bty = "n")




# To get the final epidemic size
View(x_res_int)
5364499 - mean(x_res_int[500, 2:501])
sd(x_res_int[500, 2:501])


