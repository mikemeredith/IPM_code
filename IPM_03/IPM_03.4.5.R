# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

# Run time approx. 1 min

library(IPMbook) ; library(jagsUI)

# 3.4 Analysis of matrix population models with Markov Chain
#     Monte Carlo (MCMC) software
# ==========================================================

# 3.4.5 Analysis of a matrix population model with multiple sources of
#       stochasticity and parameter uncertainty
# ---------------------------------------------------------------------

# Define mean, measurement error and temporal variability of the demographic parameters
mean.sj <- 0.3          # Mean juv. survival on probability scale
se.sj.e <- 0.005        # Uncertainty of mean sj as SE on natural scale
sd.sj.t <- 0.25         # Temporal variability of sj as SD on logit scale
mean.sa <- 0.55         # Mean ad. survival on probability scale
se.sa.e <- 0.005        # Uncertainty of mean saas SE on natural scale
sd.sa.t <- 0.07         # Temporal variability of sa as SD on logit scale
mean.f1 <- 1.3          # Mean productivity of 1y females on natural scale
se.f1.e <- 0.05         # Uncertainty of f1 as SE on natural scale
sd.f1.t <- 0.3          # Temporal variability of f1 as SD on natural scale
mean.fa <- 1.8          # Mean productivity of adult females on natural scale
se.fa.e <- 0.03         # Uncertainty of fa as SE on natural scale
sd.fa.t <- 0.3          # Temporal variability of fa as SD on natural scale

# Define the number of years with predictions and the Monte Carlo setting
T <- 200

# Bundle data
N1 <- 10
N2 <- 10
jags.data <- list(alpha.sj=getBeta2Par(mean.sj, se.sj.e)[1], beta.sj=getBeta2Par(mean.sj,
    se.sj.e)[2], alpha.sa=getBeta2Par(mean.sa, se.sa.e)[1], beta.sa=getBeta2Par(mean.sa,
    se.sa.e)[2], mean.f1=mean.f1, mean.fa=mean.fa, se.f1.e=se.f1.e, se.fa.e=se.fa.e,
    sd.sj.t=sd.sj.t, sd.sa.t=sd.sa.t, sd.f1.t=sd.f1.t, sd.fa.t=sd.fa.t, T=T, N1=N1, N2=N2)

# Write JAGS model file
cat(file="model7.txt", "
model {
  # Use of RNG to account for measurement error (uncertainty)
  # in the means of the demographic rates on the natural scale
  mean.sj.e ~ dbeta(alpha.sj, beta.sj)
  mean.sa.e ~ dbeta(alpha.sa, beta.sa)
  mean.f1.e ~ dnorm(mean.f1, tau.f1.e)
  mean.fa.e ~ dnorm(mean.fa, tau.fa.e)

  # Calculate precision of the measurement error of productivity
  tau.f1.e <- pow(se.f1.e, -2)
  tau.fa.e <- pow(se.fa.e, -2)

  # Calculate precision for the temporal variability of demographic rates
  tau.logit.sj.t <- pow(sd.sj.t, -2)
  tau.logit.sa.t <- pow(sd.sa.t, -2)
  tau.f1.t <- pow(sd.f1.t, -2)
  tau.fa.t <- pow(sd.fa.t, -2)

  # Use of RNG to accomodate temporal variability of demographic rates
  # (process variability)
  for (t in 1:T){
    sj[t] <- ilogit(logit.sj[t])                    # Backtransformation from logit scale
    logit.sj[t] ~ dnorm(logit(mean.sj.e), tau.logit.sj.t)
    sa[t] <- ilogit(logit.sa[t])                    # Backtransformation from logit scale
    logit.sa[t] ~ dnorm(logit(mean.sa.e), tau.logit.sa.t)
    f1[t] ~ dnorm(mean.f1.e, tau.f1.t)T(0,)         # Truncation to positive values
    fa[t] ~ dnorm(mean.fa.e, tau.fa.t)T(0,)         # Truncation to positive values
  }

  # Model for initial state
  N[1,1] <- N1
  N[2,1] <- N2

  # Loop over time
  for (t in 1:T){
    # Population model
    N[1,t+1] ~ dpois(sj[t] * (f1[t] * N[1,t] + fa[t] * N[2,t]))
    N[2,t+1] ~ dbin(sa[t], (N[1,t] + N[2,t]))
    extinct[t] <- equals(N[1,t+1] + N[2,t+1], 0)    # Determines whether
        # population is still thriving (extinct = 0) or went extinct (extinct = 1)
  }
}
")

# Parameters monitored
parameters <- c("N", "extinct")

# MCMC settings
ni <- 50000; nt <- 1; nb <- 0; nc <- 1; na <- 0

# Call JAGS (ART <1 min) and summarize results
out7 <- jags(jags.data, NULL, parameters, "model7.txt", n.adapt=na, n.chains=nc, n.thin=nt,
    n.iter=ni, n.burnin=nb, DIC=FALSE)
print(out7, 4)

#                    mean          sd 2.5%   50%      97.5% overlap0 f
# N[1,1]          10.0000      0.0000   10  10.0     10.000    FALSE 1
# N[2,1]          10.0000      0.0000   10  10.0     10.000    FALSE 1
# N[1,2]           9.3800      3.7034    3   9.0     17.000    FALSE 1
# N[2,2]          10.9935      2.2543    7  11.0     15.000    FALSE 1
# [ ... output truncated ... ]
# extinct[1]       0.0000      0.0000    0   0.0      0.000    FALSE 1
# extinct[2]       0.0000      0.0000    0   0.0      0.000    FALSE 1
# extinct[3]       0.0000      0.0000    0   0.0      0.000    FALSE 1
# [ ... output truncated ... ]
# extinct[198]     0.3838      0.4863    0   0.0      1.000     TRUE 1
# extinct[199]     0.3842      0.4864    0   0.0      1.000     TRUE 1
# extinct[200]     0.3846      0.4865    0   0.0      1.000     TRUE 1

# Calculation of the stochastic population growth rate outside JAGS
dimensions <- dim(out7$sims.list$extinct)
r.annual <- matrix(NA, nrow=dimensions[2], ncol=dimensions[1])
r <- lambda <- numeric(dimensions[1])
for (s in 1:dimensions[1]){
  for (t in 1:dimensions[2]){
    # Calculate annual growth rate on log scale
    if (out7$sims.list$extinct[s,t] == 1) break
    r.annual[t,s] <- log(out7$sims.list$N[s,1,t+1] + out7$sims.list$N[s,2,t+1]) -
    log(out7$sims.list$N[s,1,t] + out7$sims.list$N[s,2,t])
  } #t
  r[s] <- mean(r.annual[which(out7$sims.list$extinct[s,] == 0),s])
  lambda[s] <- exp(r[s])
} #s

mean(r)
# [1] -0.01570117
sd(r)
# [1] 0.06555304

# Population growth rate of populations that did not go extinct
mean(r[out7$sims.list$extinct[,T]==0])
# [1] 0.02353183
sd(r[out7$sims.list$extinct[,T]==0])
# [1] 0.01234406
