# Trying to construct the SIR model and simulate data from it (with intervention)

# Importing a package for differential equations - can be extended to discrete time as well
library(odin)
library(dde)
library(scales)

# i indicates which simulation
# Generating the Discrete Stochastic SEIR model without intervention dynamics
sir_mult_generator <- odin::odin({
  ## Core equations for transitions between compartments:
  update(S[]) <- S[i] - C[i]
  update(I[]) <- I[i] + C[i] - D[i]
  update(R[]) <- R[i] + D[i]
  
  ## Individual probabilities of transition: (intervention part added here)
  p_SI[] <- if (step>130) 1 - exp(-beta*exp(-q*(step-130))*I[i]/N[i]) else 1 - exp(-beta*I[i]/N[i])  # S to I
  p_IR <- 1 - exp(-gamma)     # I to R
  
  ## Draws from binomial distributions for numbers changing between compartments:
  C[] <- rbinom(S[i], p_SI[i])
  D[] <- rbinom(I[i], p_IR)
  
  ## Total population size
  N[] <- S[i] + I[i] + R[i]
  
  ## Initial states:
  initial(S[]) <- S_ini
  initial(I[]) <- I_ini
  initial(R[]) <- R_ini
  
  ## User defined parameters - default in parentheses:
  S_ini <- user(5364499)
  I_ini <- user(1)
  R_ini <- user(0)
  beta <- user(0.10008035)
  gamma <- user(0.09375407)
  q <- user(0.07254123)
  
  ## Number of replicates
  nsim <- user(500)
  dim(N) <- nsim
  dim(S) <- nsim
  dim(I) <- nsim
  dim(R) <- nsim
  dim(p_SI) <- nsim
  dim(C) <- nsim
  dim(D) <- nsim
  
}, verbose = FALSE)


# Executing the simulation of a scenario


# Running the model simulation - if you want different parameters, put them here, otherwise they will 
# be taken from inside the function above, as defined by the user
sir <- sir_mult_generator()
sir
x <- sir_mult_generator()


# Generating a graph


set.seed(50)
sir_col <- c("gold", "maroon", "sky blue")
#sir_col_transp <- paste0(sir_col, "66")
x_res_int <- x$run(0:500)
par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
matplot(x_res_int[, 1], x_res_int[, -1][,1:1000], xlab = "Time", ylab = "Number of individuals",
        type = "l", col = alpha(rep(sir_col, each = 500), 0.5), lty = 1, mgp=c(3.2,0.4,0), 
        xlim=c(0,300), ylim = c(0,130))
legend("topright", lwd = 1, col = sir_col, legend = c("S", "I", "R"), bty = "n")




# To get the final epidemic size
View(x_res_int)

# Now, we need to get the means of S,I,R and SD - confidence intervals


# Based on susceptable definition - s0 - s(end time)
5364499 - mean(x_res_int[500, 2:501])

sd(x_res_int[500, 2:501])
