# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

# 3.3 Classical analysis of a matrix population model
# ===================================================

# 3.3.3 Analysis of a matrix population model with environmental stochasticity
# ----------------------------------------------------------------------------

# Define mean and temporal variability (SD) of the demographic parameters
mean.sj <- 0.3          # Mean juvenile survival (probability scale)
sd.sj.t <- 0.25         # Temporal variability on the logit scale
mean.sa <- 0.55         # Mean adult survival (probability scale)
sd.sa.t <- 0.07         # Temporal variability on the logit scale
mean.f1 <- 1.3          # Mean productivity of 1y old females
sd.f1.t <- 0.3          # Temporal variability on the natural scale
mean.fa <- 1.8          # Mean productivity of adult females
sd.fa.t <- 0.3          # Temporal variability on the natural scale

# Define the number of years with predictions and burn-in length
T <- 100000             # Length of Markov chain
u <- 1000               # Length of burn-in period

# Generate demographic values from normal distributions
sj <- plogis(rnorm(T, qlogis(mean.sj), sd.sj.t))
sa <- plogis(rnorm(T, qlogis(mean.sa), sd.sa.t))
f1 <- rnorm(T, mean.f1, sd.f1.t)
fa <- rnorm(T, mean.fa, sd.fa.t)

# Define population matrix and initial stage-specific population sizes
N <- matrix(NA, nrow=2, ncol=T+1)
N[,1] <- c(1, 1)

# Project population forwards
r <- numeric(T)
for (t in 1:T){
  if(t %% 1000 == 0) {cat(paste("*** Simrep", t, "***\n")) }  # Counter
  A <- matrix(c(sj[t] * f1[t], sj[t] * fa[t], sa[t], sa[t]), ncol=2, byrow=TRUE)
  N[,t+1] <- A %*% N[,t]
  r[t] <- log(sum(N[,t+1])) - log(sum(N[,t]))                 # Annual population growth rate
  N[,t+1] <- N[,t+1] / sum(N[,t+1])                           # Scale N to avoid numerical overflow
}

mean(r[u:T])
# [1] 0.01944997

exp(mean(r[u:T]))
# [1] 1.01964

# Generate demographic values from normal distributions
sj <- plogis(rnorm(T, qlogis(mean.sj), sd.sj.t))
sa <- plogis(rnorm(T, qlogis(mean.sa), sd.sa.t))
f1 <- rnorm(T, mean.f1, sd.f1.t)
fa <- rnorm(T, mean.fa, sd.fa.t)

# Define population matrix and initial stage-specific population sizes
N <- N.star <- matrix(NA, nrow=2, ncol=T+1)
N[,1] <- N.star[,1] <- c(1, 1)

# Project population and calculate stochastic population growth rate
delta <- 0.001                                                # Magnitude of perturbation (= "small change")
r <- r.star <- numeric(T)
for (t in 1:T){
  if(t %% 1000 == 0) {cat(paste("*** Simrep", t, "***\n")) }  # Counter
  # Projection using sj
  A <- matrix(c(sj[t] * f1[t], sj[t] * fa[t], sa[t], sa[t]), ncol=2, byrow=TRUE)
  N[,t+1] <- A %*% N[,t]
  r[t] <- log(sum(N[,t+1])) - log(sum(N[,t]))                 # Annual population growth rate
  N[,t+1] <- N[,t+1] / sum(N[,t+1])                           # Scale N to avoid numerical overflow
  # Projection using sj.star
  A.star <- matrix(c((sj[t] + delta) * f1[t], (sj[t] + delta) * fa[t], sa[t], sa[t]), ncol=2,
      byrow=TRUE)
  N.star[,t+1] <- A.star %*% N.star[,t]
  r.star[t] <- log(sum(N.star[,t+1])) - log(sum(N.star[,t]))  # Annual population growth rate
  N.star[,t+1] <- N.star[,t+1] / sum(N.star[,t+1])            # Scale N.star to avoid numerical overflow
}

# Compute stochastic sensitivity (juvenile survival)
(exp(mean(r.star[u:T])) - exp(mean(r[u:T]))) / delta
# [1] 1.450347

# Compute stochastic elasticity (juvenile survival)
(exp(mean(r.star[u:T])) - exp(mean(r[u:T]))) / delta * mean.sj / exp(mean(r[u:T]))
# [1] 0.4263838
