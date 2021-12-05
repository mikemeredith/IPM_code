# Schaub & Kéry (2022) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------

library(IPMbook) ; library(jagsUI)

# ~~~ this needs data created in section 4.4.1 ~~~
nbrood <- 1000           # Number of broods with young counted
brood.mean <- 1.5        # Average brood size
sd.brood <- 0.3          # log-linear brood random effect
set.seed(24)
expNyoung <- exp(log(brood.mean) + rnorm(nbrood, 0, sd.brood))
C <- rpois(nbrood, expNyoung)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 4.4 Models for productivity surveys
# ===================================

# 4.4.4 Censoring in brood size data
# ----------------------------------

# Simulate right-censored Poisson observations
set.seed(1)
C3 <- C                                         # Make another copy of the data set
sum(C3 >= 2)                                    # How many broods 2 or greater ? -- 439
large.broods <- which(C3>=2)                    # Index them
censored.broods <- sample(large.broods, 300)    # Randomly select 300
C3[censored.broods] <- C3[censored.broods] - 1  # Censor these 300
mean(C)                                         # True mean
mean(C3)                                        # Mean of C3 biased low
# [1] 1.529
# [1] 1.229

C4 <- C3                                        # Make another copy
C4[censored.broods] <- NA                       # Censored are NA'd out
d <- as.numeric(is.na(C4))                      # Binary censoring indicator

# Compute censoring threshold
# Equal to observed value for censored
# Equal to observed value plus some small value for non-censored
cens.threshold <- C3
cens.threshold[d==0] <- C3[d==0] + 0.1
cbind('Truth'=C, 'Obs w/cens'=C3, 'Modelled data (C4)'=C4, 'Cens. indicator (d)'=d, cens.threshold)[1:6,]
              # Look at first 6 rows
#      Truth Obs w/cens Modelled data (C4) Cens. indicator (d) cens.threshold
# [1,]     0          0                  0                   0            0.1
# [2,]     0          0                  0                   0            0.1
# [3,]     1          1                  1                   0            1.1
# [4,]     1          1                  1                   0            1.1
# [5,]     4          3                 NA                   1            3.0
# [6,]     1          1                  1                   0            1.1

# Data bundle
jags.data <- list(C4=C4, d=d, cens.threshold=cens.threshold, nbrood=nbrood)
str(jags.data)
# List of 4
# $ C4            : num [1:1000] 0 0 1 1 NA 1 0 NA 0 1 ...
# $ d             : num [1:1000] 0 0 0 0 1 0 0 1 0 0 ...
# $ cens.threshold: num [1:1000] 0.1 0.1 1.1 1.1 3 1.1 0.1 1 0.1 1.1 ...
# $ nbrood        : num 1000

# Write JAGS model file
cat(file="model11.txt", "
model {
  # Prior
  rho ~ dunif(0, 5)

  # Likelihood
  for (i in 1:nbrood){
    C4[i] ~ dpois(rho)
    prob.censored[i] <- step(C4[i] - cens.threshold[i])
    d[i] ~ dbern(prob.censored[i])
  }
}
")

# Initial values
C4st <- rep(NA, nbrood)
C4st[is.na(C4)] <- C3[is.na(C4)]
inits <- function() list(rho=runif(1), C4=C4st)

# Parameters monitored
parameters <- c("rho", "C4", "prob.censored")

# MCMC settings
ni <- 2000; nb <- 1000; nc <- 3; nt <- 2; na <- 1000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out14 <- jags(jags.data, inits, parameters, "model11.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out14) # Not shown
print(out14, 3)

#                         mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# rho                    1.446  0.042    1.363    1.445    1.531    FALSE 1 1.001  1500
# C4[1]                  0.000  0.000    0.000    0.000    0.000    FALSE 1    NA     1
# C4[2]                  0.000  0.000    0.000    0.000    0.000    FALSE 1    NA     1
# C4[3]                  1.000  0.000    1.000    1.000    1.000    FALSE 1    NA     1
# C4[4]                  1.000  0.000    1.000    1.000    1.000    FALSE 1    NA     1
# C4[5]                  3.459  0.733    3.000    3.000    5.000    FALSE 1 1.000  1500
# C4[6]                  1.000  0.000    1.000    1.000    1.000    FALSE 1    NA     1
# [ .... ]

# Remind ourselves of mean observed brood size ignoring censoring
mean(C3)
# [1] 1.229

C5 <- C                                         # Make another copy of the data
# Assume last 200 broods only known to be a success or failure
C5part1 <- C5[1:800]                            # The unaltered counts
C5part2 <- as.numeric(C5[801:1000]>0)           # Success/failure indicators

# Data bundle
jags.data <- list(C5part1copy=C5part1, C5part1=C5part1, C5part2=C5part2)
str(jags.data)
# List of 3
# $ C5part1copy: int [1:800] 0 0 1 1 4 1 0 2 0 1 ...
# $ C5part1    : int [1:800] 0 0 1 1 4 1 0 2 0 1 ...
# $ C5part2    : num [1:200] 0 1 0 1 1 1 1 1 0 1 ...

# Write JAGS model file
cat(file="model12.txt", "
model {
  # Priors
  rho1 ~ dunif(0, 5)                            # For model part 1
  rho2 ~ dunif(0, 5)                            # For model part 2

  # Likelihood
  # First model: Poisson GLM for C5 (part 1)
  for (i in 1:800){
    C5part1copy[i] ~ dpois(rho1)
  }

  # Second, integrated model: Poisson GLM for C5 (part 1) plus a Bernoulli GLM for C5 (part 2) with
  # appropriate definition of the success probability as a function of the same rho
  for (i in 1:800){                             # Submodel for part 1 of the data
    C5part1[i] ~ dpois(rho2)
  }
  for (i in 1:200){                             # Submodel for part 2 of the data
    C5part2[i] ~ dbern(1 - exp(-rho2))
  }
}
")

# Initial values
inits <- function() list(rho1=runif(1), rho2=runif(1))

# Parameters monitored
parameters <- c("rho1", "rho2")

# MCMC settings
ni <- 5000; nb <- 1000; nc <- 3; nt <- 4; na <- 1000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out15 <- jags(jags.data, inits, parameters, "model12.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out15) # Not shown
print(out15, 3)

#              mean    sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# rho1        1.545 0.045    1.461    1.545    1.638    FALSE 1 1.001  1173
# rho2        1.526 0.041    1.445    1.526    1.608    FALSE 1 1.001  1053
