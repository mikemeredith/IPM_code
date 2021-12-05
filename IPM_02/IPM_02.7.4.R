# Schaub & Kéry (2022) Integrated Population Models
# Chapter 2 : Bayesian statistical modeling using JAGS
# ----------------------------------------------------

library(IPMbook) ; library(jagsUI)

# ~~~ needs this from section 2.7.3 ~~~
nsites <- 1000                          # Number of sites
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2.7 Using JAGS to fit simple statistical models from R: GLMs and GLMMs
# ======================================================================

# 2.7.4 Multinomial generalized linear models
# -------------------------------------------

# Simulate a habitat data set by drawing multinomial response with stated cell probs pi
set.seed(48)
hab.types <- c("cliff top", "forest edge", "rocky pasture", "scree")
pi <- c(0.09, 0.32, 0.17, 0.42)             # Relative probabilities, must sum to 1
habitat <- rmultinom(1, nsites, pi)
dimnames(habitat) <- list(hab.types, "Number of sites")
habitat

#               Number of sites
# cliff top                  79
# forest edge               313
# rocky pasture             181
# scree                     427

# Bundle data
jags.data <- list(C=as.numeric(habitat), k=length(habitat), N=sum(habitat))
str(jags.data)

# List of 3
# $ C: num [1:4] 79 313 181 427
# $ k: int 4
# $ N: int 1000

# Write JAGS model file
cat(file="model4.txt", "
model {
  # Priors and linear models
  # Cell probabilities for the multinomial
  for (i in 1:k){                           # Loop over all cells in the multinomial
    d[i] ~ dgamma(1, 1)
    pi[i] <- d[i] / sum(d[])
  }
  # Likelihood for set of multinomial counts
  C ~ dmulti(pi, N)
}
")

# Initial values
inits <- function(){list(d=rgamma(4, 1, 1))}

# Parameters monitored
parameters <- c("pi")

# MCMC settings
ni <- 110000; nb <- 10000; nc <- 3; nt <- 10; na <- 1000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out5 <- jags(jags.data, inits, parameters, "model4.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out5) # Not shown
print(out5, 3)

#            mean    sd   2.5%    50%  97.5% overlap0 f Rhat n.eff
# pi[1]     0.080 0.009  0.064  0.079  0.097    FALSE 1    1 30000
# pi[2]     0.313 0.015  0.284  0.313  0.342    FALSE 1    1 12397
# pi[3]     0.181 0.012  0.158  0.181  0.206    FALSE 1    1 16620
# pi[4]     0.426 0.016  0.396  0.426  0.457    FALSE 1    1 30000

# ~~~~ save output for use later ~~~~
save(out5, file="IPM_02.7.4_out5.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Simulate survival data: individual encounter histories
library(AHMbook)
set.seed(39)
dat <- simCJS(n.occ=4, n.marked=30, phi=0.7, p=0.4, show.plot=TRUE)
dat$ch                                      # Look at these encounter history data (not shown)

# Produce the m-array from encounter histories
Cmat <- marray(dat$ch)
Cmat
#         recaptured
# released Y2 Y3 Y4 never
#       Y1  7  7  1    15
#       Y2  0 12  5    20
#       Y3  0  0 13    36

# Compute size of release cohorts
R <- rowSums(Cmat)

# Bundle data
jags.data <- list(marray=Cmat, R=R)
str(jags.data)
# List of 2
# $ marray: num [1:3, 1:4] 7 0 0 7 12 0 1 5 13 15 ...
# $ R     : num [1:3] 30 37 49

# Write JAGS model file
cat(file="model5.txt", "
model {
  # Priors and linear models
  phi ~ dunif(0, 1)                         # Apparent survival probability
  p ~ dunif(0, 1)                           # Recapture probability

  # Vectors of cell probabilities in the three multinomials
  # Release cohort 1
  pi1[1] <- phi * p                         # First recaptured in recap year 2
  pi1[2] <- phi^2 * (1-p) * p               # First recaptured in recap year 3
  pi1[3] <- phi^3 * (1-p)^2 * p             # First recaptured in recap year 4
  pi1[4] <- 1 - sum(pi1[1:3])               # Never recaptured

  # Release cohort 2
  pi2[1] <- 0                               # Takes account of structural zero counts
  pi2[2] <- phi * p
  pi2[3] <- phi^2 * (1-p) * p
  pi2[4] <- 1 - sum(pi2[1:3])

  # Release cohort 3
  pi3[1] <- 0                               # Accounts for structural zero counts
  pi3[2] <- 0                               # Accounts for structural zero counts
  pi3[3] <- phi * p
  pi3[4] <- 1 - sum(pi3[1:3])

  # Likelihood: one multinomial for each row in the m-array
  marray[1,] ~ dmulti(pi1, R[1])            # MN for release cohort 1
  marray[2,] ~ dmulti(pi2, R[2])            # MN for release cohort 2
  marray[3,] ~ dmulti(pi3, R[3])            # MN for release cohort 3
}
")

# Initial values
inits <- function(){list(phi=runif(1, 0, 1), p=runif(1, 0, 1))}

# Parameters monitored
parameters <- c("phi", "p", "pi1", "pi2", "pi3")

# MCMC settings
ni <- 10000; nb <- 5000; nc <- 3; nt <- 5; na <- 1000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out6 <- jags(jags.data, inits, parameters, "model5.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out6) # Not shown
print(out6, 3)

#            mean    sd   2.5%    50%  97.5% overlap0 f  Rhat n.eff
# phi       0.786 0.095  0.603  0.783  0.967    FALSE 1 1.001  2273
# p         0.376 0.079  0.241  0.369  0.540    FALSE 1 1.000  3000
# pi1[1]    0.290 0.041  0.212  0.288  0.373    FALSE 1 1.000  3000
# pi1[2]    0.141 0.026  0.089  0.141  0.191    FALSE 1 1.000  2950
# pi1[3]    0.072 0.027  0.027  0.070  0.129    FALSE 1 1.001  2223
# pi1[4]    0.497 0.061  0.375  0.497  0.617    FALSE 1 1.001  3000
# [... output truncated ...]

# Print the pi matrix: cell probabilities for all 3 release cohorts
print(rbind(out6$mean$pi1, out6$mean$pi2, out6$mean$pi3), 3)
#      [,1]  [,2]   [,3]  [,4]
# [1,] 0.29 0.141 0.0723 0.497
# [2,] 0.00 0.290 0.1410 0.569
# [3,] 0.00 0.000 0.2896 0.710

print(mean(out6$sims.list$phi * out6$sims.list$p), 3)
# [1] 0.29

# ~~~~ save output for use later ~~~~
save(out6, file="IPM_02.7.4_out6.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code  for MLE analysis ~~~~
# For a quick comparison with MLEs of the same model, we use functionality in Mike Meredith's wiqid package and find slightly less agreement than what we have perhaps come to expect. The reasons for which may be small sample size and asymmetric posterior distributions.

library(wiqid)
?survCJS                 # Check out to see how it works
(mle <- survCJS(dat$ch)) # Fit phi(.)p(.) variant of CJS

# Real values (duplicates omitted):
#         est  lowCI  uppCI
# phi1 0.7938 0.5137 0.9334
# p1   0.3610 0.2188 0.5325
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
