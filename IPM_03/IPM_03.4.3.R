# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

library(IPMbook) ; library(jagsUI)

# 3.4 Analysis of matrix population models with Markov Chain
#     Monte Carlo (MCMC) software
# ==========================================================

# 3.4.3 Analysis of a matrix population model with environmental stochasticity
# ----------------------------------------------------------------------------

# Define mean and temporal variability (SD) of the demographic parameters
mean.sj <- 0.3            # Mean juv. survival on probability scale
sd.sj.t <- 0.25           # Temporal variability of sj on the logit scale
mean.sa <- 0.55           # Mean ad. survival on probability scale
sd.sa.t <- 0.07           # Temporal variability of sa on the logit scale
mean.f1 <- 1.3            # Mean productivity of 1y females on natural scale
sd.f1.t <- 0.3            # Temporal variability of f1 on the natural scale
mean.fa <- 1.8            # Mean productivity of adult females on natural scale
sd.fa.t <- 0.3            # Temporal variability of fa on the natural scale

# Bundle data
jags.data <- list(mean.sj=mean.sj, mean.sa=mean.sa, mean.f1=mean.f1, mean.fa=mean.fa,
    sd.sj.t=sd.sj.t, sd.sa.t=sd.sa.t, sd.f1.t=sd.f1.t, sd.fa.t=sd.fa.t, T=21000, u=1000)

# Write JAGS model file
cat(file="model4.txt", "
model {
  # Calculate precision for temporal variability of demographic rates
  tau.logit.sj <- pow(sd.sj.t, -2)
  tau.logit.sa <- pow(sd.sa.t, -2)
  tau.f1 <- pow(sd.f1.t, -2)
  tau.fa <- pow(sd.fa.t, -2)

  # Use of RNG to accomodate temporal variability of demographic rates (process variability)
  for (t in 1:T){
    sj[t] <- ilogit(logit.sj[t])                      # Backt. from logit to natural scale
    logit.sj[t] ~ dnorm(logit(mean.sj), tau.logit.sj)
    sa[t] <- ilogit(logit.sa[t])                      # Backt. from logit to natural scale
    logit.sa[t] ~ dnorm(logit(mean.sa), tau.logit.sa)
    f1[t] ~ dnorm(mean.f1, tau.f1)
    fa[t] ~ dnorm(mean.fa, tau.fa)
  }

  # Model for initial state
  N[1,1] <- 1
  N[2,1] <- 1

  # Loop over time
  for (t in 1:T){
    # Population model
    N[1,t+1] <- sj[t] * (f1[t] * N[1,t] + fa[t] * N[2,t])
    N[2,t+1] <- sa[t] * (N[1,t] + N[2,t])

    # Annual growth rate on log scale
    r.annual[t] <- log(N[1,t+1] + N[2,t+1]) - log(N[1,t] + N[2,t])
  }
  r <- mean(r.annual[u:T])
  lambda <- exp(r)

  # Sensitivity and elasticity of lambda to changes in sj
  delta <- 0.001                                      # Size of perturbation
  N.star[1,1] <- 1
  N.star[2,1] <- 1
  for (t in 1:T){
    N.star[1,t+1] <- (sj[t] + delta) * (f1[t] * N.star[1,t] + fa[t] * N.star[2,t])
    N.star[2,t+1] <- sa[t] * (N.star[1,t] + N.star[2,t])
    r.annual.star[t] <- log(N.star[1,t+1] + N.star[2,t+1]) - log(N.star[1,t] + N.star[2,t])
  }
  r.star <- mean(r.annual.star[u:T])
  s.sj <- (exp(r.star) - lambda) / delta
  e.sj <- s.sj * mean.sj / lambda
}
")

# Parameters monitored
parameters <- c("r", "lambda", "s.sj", "e.sj")

# MCMC settings
ni <- 1; nt <- 1; nb <- 0; nc <- 1; na <- 0

# Call JAGS (ART <1 min) and summarize results
out4 <- jags(jags.data, NULL, parameters, "model4.txt", n.adapt=na, n.chains=nc, n.thin=nt, n.iter=ni,
    n.burnin=nb, DIC=FALSE)
print(out4, 4)

#          mean sd   2.5%    50%  97.5% overlap0 f
# r      0.0191 NA 0.0191 0.0191 0.0191    FALSE 1
# lambda 1.0193 NA 1.0193 1.0193 1.0193    FALSE 1
# s.sj   1.4486 NA 1.4486 1.4486 1.4486    FALSE 1
# e.sj   0.4263 NA 0.4263 0.4263 0.4263    FALSE 1


sj.t <- c(0.303, 0.300, 0.288, 0.310, 0.306, 0.291, 0.303, 0.281, 0.314, 0.303)
sa.t <- c(0.546, 0.545, 0.547, 0.545, 0.547, 0.556, 0.540, 0.551, 0.548, 0.550)
f1.t <- c(1.56, 1.15, 1.38, 1.89, 1.06, 0.69, 1.16, 1.10, 1.48, 1.28)
fa.t <- c(2.03, 1.69, 1.95, 1.98, 1.63, 1.43, 1.63, 1.45, 1.79, 1.91)

T <- 21000                                            # Number of predicted years
t <- 1:10                                             # Number of years with actual data
year <- sample(t, T, replace=TRUE)
sj <- sj.t[year]
sa <- sa.t[year]
f1 <- f1.t[year]
fa <- fa.t[year]

# Bundle data
jags.data <- list(sj=sj, sa=sa, f1=f1, fa=fa, T=21000, u=1000)

# Write JAGS model file
cat(file="model5.txt", "
model {
  # Model for initial state
  N[1,1] <- 1
  N[2,1] <- 1

  # Loop over time
  for (t in 1:T){
    # Population model
    N[1,t+1] <- sj[t] * (f1[t] * N[1,t] + fa[t] * N[2,t])
    N[2,t+1] <- sa[t] * (N[1,t] + N[2,t])

    # Annual (realized) growth rate on log scale
    r.annual[t] <- log(N[1,t+1] + N[2,t+1]) - log(N[1,t] + N[2,t])
  }
  r <- mean(r.annual[u:T])
  lambda <- exp(r)

  # Sensitivity and elasticity of lambda to changes in sj
  delta <- 0.001                                                # size of perturbation
  N.star[1,1] <- 1
  N.star[2,1] <- 1
  for (t in 1:T){
    N.star[1,t+1] <- (sj[t] + delta) * (f1[t] * N.star[1,t] + fa[t] * N.star[2,t])
    N.star[2,t+1] <- sa[t] * (N.star[1,t] + N.star[2,t])
    r.annual.star[t] <- log(N.star[1,t+1] + N.star[2,t+1]) - log(N.star[1,t] + N.star[2,t])
  }
  r.star <- mean(r.annual.star[u:T])
  s.sj <- (exp(r.star) - lambda) / delta
  e.sj <- s.sj * mean(sj) / lambda
}
")

# Parameters monitored
parameters <- c("r", "lambda", "s.sj", "e.sj")

# MCMC settings
ni <- 1; nt <- 1; nb <- 0; nc <- 1; na <- 0

# Call JAGS (ART <1 min) and summarize results
out5 <- jags(jags.data, NULL, parameters, "model5.txt", n.adapt=na, n.chains=nc, n.thin=nt,
    n.iter=ni, n.burnin=nb, DIC=FALSE)
print(out5, 4)

#          mean sd   2.5%    50%  97.5% overlap0 f
# r      0.0061 NA 0.0061 0.0061 0.0061    FALSE 1
# lambda 1.0061 NA 1.0061 1.0061 1.0061    FALSE 1
# s.sj   1.4112 NA 1.4112 1.4112 1.4112    FALSE 1
# e.sj   0.4207 NA 0.4207 0.4207 0.4207    FALSE 1
