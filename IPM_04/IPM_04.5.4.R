# Schaub & Kéry (2022) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------

# Run time approx. 80 secs

# 4.5 Models for survival surveys
# ===============================

# 4.5.4 Joint analysis of capture-recapture and dead-recovery data
# ----------------------------------------------------------------

library(IPMbook); library(jagsUI)
data(stork)
str(stork)
# int [1:691, 1:16] 0 0 0 0 0 0 0 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : NULL
# ..$ : chr [1:16] "1986" "1987" "1988" "1989" ...

marr <- marray(stork, unobs=2)

# Bundle data
unobs <- 2                                            # Number of unobserved states
ns <- length(unique(as.numeric(stork))) - 1 + unobs   # Total num. states
jags.data <- list(marr=marr, nyears=ncol(stork), rel=rowSums(marr), ns=ns,
    zero=matrix(0, ns, ns), ones = diag(ns))
str(jags.data)
# List of 6
# $ marr  : num [1:60, 1:61] 1 0 0 0 0 0 0 0 0 0 ...
# $ nyears: int 16
# $ rel   : num [1:60] 4 0 0 0 15 0 0 0 30 0 ...
# $ ns    : num 4
# $ zero  : num [1:4, 1:4] 0 0 0 0 0 0 0 0 0 0 ...
# $ ones  : num [1:4, 1:4] 1 0 0 0 0 1 0 0 0 0 ...

# Write JAGS model file
cat(file="model22.txt", "
model {
  # Priors and linear models
  for (t in 1:(nyears-1)){
    logit.s[t] ~ dnorm(mu, tau)                       # Means parameterization of random effect
    s[t] <- ilogit(logit.s[t])
    f[t] <- mean.f
    r[t] <- mean.r
    p[t] ~ dunif(0, 1)
  }

  mean.s ~ dunif(0, 1)                                # Prior for mean survival
  mean.f ~ dunif(0, 1)                                # Prior for mean site fidelity
  mean.r ~ dunif(0, 1)                                # Prior for mean recovery
  mu <- logit(mean.s)
  sigma ~ dunif(0, 10)                                # Prior for temporal variability of survival
  tau <- pow(sigma, -2)

  # Define state-transition and observation probabilities
  for (t in 1:(nyears-1)){
    psi[1,t,1] <- s[t] * f[t]                         # State-transitions
    psi[1,t,2] <- 1-s[t]
    psi[1,t,3] <- s[t] * (1-f[t])
    psi[1,t,4] <- 0
    psi[2,t,1] <- 0
    psi[2,t,2] <- 0
    psi[2,t,3] <- 0
    psi[2,t,4] <- 1
    psi[3,t,1] <- 0
    psi[3,t,2] <- 1-s[t]
    psi[3,t,3] <- s[t]
    psi[3,t,4] <- 0
    psi[4,t,1] <- 0
    psi[4,t,2] <- 0
    psi[4,t,3] <- 0
    psi[4,t,4] <- 1

    po[1,t] <- p[t]                                   # Observation probabilities
    po[2,t] <- r[t]
    po[3,t] <- 0
    po[4,t] <- 0
  }

  # From here onwards, no changes needed regardless of which model is fitted
  # Calculate probability of non-encounter (dq) and reshape the array for the encounter
  # probabilities
  for (t in 1:(nyears-1)){
    for (s in 1:ns){
      dp[s,t,s] <- po[s,t]
      dq[s,t,s] <- 1-po[s,t]
    } #s
    for (s in 1:(ns-1)){
      for (m in (s+1):ns){
        dp[s,t,m] <- 0
        dq[s,t,m] <- 0
      } #s
    } #m
    for (s in 2:ns){
      for (m in 1:(s-1)){
        dp[s,t,m] <- 0
        dq[s,t,m] <- 0
      } #s
    } #m
  } #t

  # Define the multinomial likelihood
  for (t in 1:((nyears-1)*ns)){
    marr[t,1:(nyears *ns-(ns-1))] ~ dmulti(pi[t,], rel[t])
  }

  # Define cell probabilities of the multistate m-array
  # Matrix U: product of probabilities of state-transition and non-encounter (needed because
  # there is no product function for matrix multiplication in JAGS)
  for (t in 1:(nyears-2)){
    U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- ones
    for (j in (t+1):(nyears-1)){
      U[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(j-2)*ns+(1:ns)] %*% psi[,t,] %*%
          dq[,t,]
    } #j
  } #t
  U[(nyears-2)*ns+(1:ns), (nyears-2)*ns+(1:ns)] <- ones

  for (t in 1:(nyears-2)){
    # Diagonal
    pi[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% psi[,t,] %*%
        dp[,t,]
    # Above main diagonal
    for (j in (t+1):(nyears-1)){
      pi[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] %*% psi[,j,] %*%
          dp[,j,]
    } #j
  } #t
  pi[(nyears-2)*ns+(1:ns),(nyears-2)*ns+(1:ns)] <- psi[,nyears-1,] %*% dp[,nyears-1,]

  # Below main diagonal
  for (t in 2:(nyears-1)){
    for (j in 1:(t-1)){
      pi[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
    } #j
  } #t

  # Last column: probability of non-recapture
  for (t in 1:((nyears-1)*ns)){
    pi[t,(nyears*ns-(ns-1))] <- 1-sum(pi[t,1:((nyears-1)*ns)])
  } #t
}
")

# Initial values
inits <- function(){list(mean.f=runif(1, 0, 1))}

# Parameters monitored
parameters <- c("mean.s", "mean.f", "mean.r", "sigma", "s", "p")

# MCMC settings
ni <- 10000; nb <- 5000; nc <- 3; nt <- 1; na <- 5000

# Call JAGS from R (ART 1 min) and check convergence
out25 <- jags(jags.data, inits, parameters, "model22.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
# par(mfrow=c(3, 3)); traceplot(out25) # Not shown
traceplot(out25) # Not shown
print(out25)
#             mean    sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# mean.s     0.793 0.047   0.717   0.788   0.890    FALSE 1 1.004   716
# mean.f     0.907 0.049   0.809   0.911   0.989    FALSE 1 1.004   660
# mean.r     0.043 0.011   0.025   0.041   0.068    FALSE 1 1.007  1047
# sigma      0.291 0.158   0.021   0.272   0.667    FALSE 1 1.015   180
# s[1]       0.792 0.066   0.660   0.791   0.916    FALSE 1 1.002  1217
# [... output truncated ...]
# s[15]      0.798 0.065   0.669   0.797   0.924    FALSE 1 1.001  3594
# p[1]       0.488 0.230   0.112   0.461   0.945    FALSE 1 1.000 10916
# [... output truncated ...]
# p[15]      0.720 0.071   0.593   0.716   0.870    FALSE 1 1.002  2696
