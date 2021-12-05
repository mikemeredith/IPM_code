# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

# Run time approx. 1 min

library(IPMbook) ; library(jagsUI)

# 3.4 Analysis of matrix population models with Markov Chain
#     Monte Carlo (MCMC) software
# ==========================================================

# 3.4.4 Analysis of a matrix population model with demographic stochasticity
# --------------------------------------------------------------------------

# Define mean of the demographic parameters
mean.sj <- 0.3
mean.sa <- 0.55
mean.f1 <- 1.3
mean.fa <- 1.8

# Define the number of years with predictions and the Monte Carlo setting
T <- 200

# Define initial stage-specific population sizes
N1 <- 10
N2 <- 10

# Bundle data
jags.data <- list(sj=mean.sj, sa=mean.sa, f1=mean.f1, fa=mean.fa, T=T, N1=N1, N2=N2)

# Write JAGS model file
cat(file="model6.txt", "
model {
  # Model for initial state
  N[1,1] <- N1
  N[2,1] <- N2

  # Loop over time
  for (t in 1:T){
    # Population model
    N[1,t+1] ~ dpois(sj * (f1 * N[1,t] + fa * N[2,t]))
    N[2,t+1] ~ dbin(sa, (N[1,t] + N[2,t]))
    extinct[t] <- equals(N[1,t+1] + N[2,t+1], 0)      # Determines whether
        # population is still thriving (extinct = 0) or went extinct (extinct = 1)
  }
}
")

# Initial values
Ninit <- matrix(10, nrow=2, ncol=T+1)
Ninit[,1] <- NA
inits <- function(){list(N=Ninit)}

# Parameters monitored
parameters <- c("N", "extinct")

# MCMC settings
ni <- 50000; nt <- 1; nb <- 0; nc <- 1; na <- 0

# Call JAGS (ART <1 min) and summarize results
out6 <- jags(jags.data, inits, parameters, "model6.txt", n.chains=nc, n.thin=nt, n.iter=ni,
    n.burnin=nb, DIC=FALSE)
print(out6, 4)

#                  mean       sd 2.5%   50%    97.5% overlap0 f
# N[1,1]        10.0000   0.0000   10  10.0   10.000    FALSE 1
# N[2,1]        10.0000   0.0000   10  10.0   10.000    FALSE 1
# N[1,2]         9.3110   3.0491    4   9.0   16.000    FALSE 1
# N[2,2]        10.9922   2.2313    7  11.0   15.000    FALSE 1
# [ ...output truncated... ]
# extinct[1]     0.0000   0.0000    0   0.0    0.000    FALSE 1
# extinct[2]     0.0000   0.0000    0   0.0    0.000    FALSE 1
# extinct[3]     0.0000   0.0000    0   0.0    0.000    FALSE 1
# [output truncated... ]
# extinct[198]   0.2876   0.4526    0   0.0    1.000     TRUE 1
# extinct[199]   0.2876   0.4527    0   0.0    1.000     TRUE 1
# extinct[200]   0.2878   0.4527    0   0.0    1.000     TRUE 1

# Calculation of the stochastic population growth rate and of the extinction probability
# outside JAGS
dimensions <- dim(out6$sims.list$extinct)
r.annual <- matrix(NA, nrow=dimensions[2], ncol=dimensions[1])
r <- lambda <- numeric(dimensions[1])
for (s in 1:dimensions[1]){
  for (t in 1:dimensions[2]){
    # Calculate annual growth rate on log scale
    if (out6$sims.list$extinct[s,t] == 1) break
    r.annual[t,s] <- log(out6$sims.list$N[s,1,t+1] + out6$sims.list$N[s,2,t+1]) -
    log(out6$sims.list$N[s,1,t] + out6$sims.list$N[s,2,t])
  } #t
  r[s] <- mean(r.annual[which(out6$sims.list$extinct[s,] == 0),s])
  lambda[s] <- exp(r[s])
} #s

mean(r)
# [1] -0.008533999
sd(r)
# [1] 0.05592566

mean(r[out6$sims.list$extinct[,T]==0])
# [1] 0.01973898
sd(r[out6$sims.list$extinct[,T]==0])
# [1] 0.006009628
