# Schaub & Kéry (2022) Integrated Population Models
# Chapter 2 : Bayesian statistical modeling using JAGS
# ----------------------------------------------------

# Run time approx. 5 mins

library(IPMbook) ; library(jagsUI)

# ~~~ requires data created in 2.7.4 ~~~
nsites <- 1000                          # Number of sites
set.seed(48)
hab.types <- c("cliff top", "forest edge", "rocky pasture", "scree")
pi <- c(0.09, 0.32, 0.17, 0.42)         # Relative probabilities, sum to 1
habitat <- rmultinom(1, nsites, pi)
dimnames(habitat) <- list(hab.types, "Number of sites")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 2.7 Using JAGS to fit simple statistical models from R: GLMs and GLMMs
# ======================================================================

# 2.7.5 Categorical generalized linear models
# -------------------------------------------

# 'Stretch out' the habitat data from Section 2.7.4.
habitat                                     # Check out habitat again (not shown)
habitat2 <- rep(1:4, habitat)               # Disaggregated version of the habitat data
habitat2; table(habitat2)                   # Look at data (not shown)

# Bundle data
jags.data <- list(C=habitat2, k=length(unique(habitat2)), n=length(habitat2))
str(jags.data)
# List of 3
# $ C: int [1:1000] 1 1 1 1 1 1 1 1 1 1 ...
# $ k: int 4
# $ n: int 1000

# Write JAGS model file
cat(file="model6.txt", "
model {
  # Priors and linear models
  # Cell probabilities
  for (i in 1:k){                           # Loop over all cells
    d[i] ~ dgamma(1, 1)
    pi[i] <- d[i] / sum(d[])
  }
  # Likelihood
  for (i in 1:n){
    C[i] ~ dcat(pi)
  }
}
")

# Initial values
inits <- function(){list(d=rgamma(4, 1, 1))}

# Parameters monitored
parameters <- c("pi")

# MCMC settings
ni <- 110000; nb <- 10000; nc <- 3; nt <- 10; na <- 1000

# Call JAGS from R (ART 5 min), check convergence and summarize posteriors
out7 <- jags(jags.data, inits, parameters, "model6.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out7) # Not shown
print(out7, 3)

#              mean    sd     2.5%      50%    97.5% overlap0 f Rhat n.eff
# pi[1]       0.080 0.009    0.064    0.079    0.097    FALSE 1    1 22081
# pi[2]       0.313 0.015    0.284    0.313    0.342    FALSE 1    1 30000
# pi[3]       0.181 0.012    0.158    0.181    0.206    FALSE 1    1 30000
# pi[4]       0.426 0.016    0.396    0.426    0.457    FALSE 1    1 14368

# Compare run times of multinomial and categorical versions of analysis
# ~~~ May need to reload out5 ~~~
load("IPM_02.7.4_out5.RData")
out7$mcmc.info$elapsed.mins / out5$mcmc.info$elapsed.mins # [1] 138.2222
