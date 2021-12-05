# Schaub & Kéry (2022) Integrated Population Models
# Chapter 2 : Bayesian statistical modeling using JAGS
# ----------------------------------------------------

library(IPMbook) ; library(jagsUI)


# 2.8 Fitting general integrated models (IMs) in JAGS
# ===================================================

# Simulate the two data sets and plot them
# Choose sample size and parameter values for both data sets
nsites1 <- 200                                      # Sample size for count data
nsites2 <- 500                                      # Sample size for detection/nondetection data
mean.lam <- 2                                       # Average expected abundance (lambda) per site
beta <- -2                                          # Coefficient of elevation covariate on lambda

# Simulate elevation covariate for both and standardize to mean of 1000 and
# standard deviation also of 1000 m
set.seed(2016)
elev1 <- sort(runif(nsites1, 200, 2000))            # Imagine 200-2000 m a.s.l.
elev2 <- sort(runif(nsites2, 200, 2000))
selev1 <- (elev1 - 1000) / 1000
selev2 <- (elev2 - 1000) / 1000

# Create two data sets with replicated counts
C1 <- rpois(nsites1, exp(log(mean.lam) + beta * selev1))
C2 <- rpois(nsites2, exp(log(mean.lam) + beta * selev2))

# Turn count data set 2 (C2) into detection/nondetection data (y)
y <- C2                                             # Make a copy
y[y > 1] <- 1                                       # Squash to binary

# ~~~~ extra code for Fig 2.19 ~~~~
# Plot counts in both data sets and detection-nondetection data
library(scales)
op <- par(mfrow=c(1, 2), mar=c(5, 5, 4, 1), cex=1.2, cex.lab=1.5, cex.axis=1.5, las=1)
plot(elev2, jitter(C2), pch=16, xlab='Elevation (m)', ylab='Counts',
    frame=FALSE, ylim = range(c(C1, C2)), col=alpha('grey80', 1))
points(elev1, jitter(C1), pch=16)
lines(200:2000, exp(log(2) -2 * ((200:2000)-1000)/1000 ), col='red', lty=1, lwd=2)
axis(1, at=c(250, 750, 1250, 1750), tcl=-0.25, labels=NA)
plot(elev2, y, xlab='Elevation (m)', ylab='Detection-Nondetection',
    axes=FALSE, pch=16, col=alpha('grey60', 0.3))
axis(1)
axis(1, at=c(250, 750, 1250, 1750), tcl=-0.25, labels=NA)
axis(2, at=c(0, 1), labels=c(0, 1))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get MLEs for individual data sets
# Data Set 1: Poisson GLM with log link for counts
summary(fm1 <- glm(C1 ~ selev1, family=poisson(link="log")))
exp(coef(fm1)[1])                                   # Estimate of lambda on natural scale from counts

# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.62035    0.06213   9.985   <2e-16 ***
# selev1      -2.19121    0.11691 -18.742   <2e-16 ***
# [  .... truncated  .....]

# (Intercept)
    # 1.85958

# Data Set 2: Bernoulli GLM with cloglog link for detection-nondetection
summary(fm2 <- glm(y ~ selev2, family=binomial(link="cloglog")))
exp(coef(fm2)[1]) # Estimate of lambda on natural scale from binary data

# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.59805    0.07968   7.506  6.1e-14 ***
# selev2      -1.95305    0.17628 -11.079  < 2e-16 ***
# [  .... truncated  .....]

# (Intercept)
   # 1.818574

# Bundle data
jags.data <- list(C=C1, y=y, nsites1=nsites1, nsites2=nsites2, selev1=selev1, selev2=selev2)
str(jags.data)
# List of 6
# $ C      : int [1:200] 14 5 11 9 8 15 11 9 9 6 ...
# $ y      : num [1:500] 1 1 1 1 1 1 1 1 1 1 ...
# $ nsites1: num 200
# $ nsites2: num 500
# $ selev1 : num [1:200] -0.795 -0.789 -0.78 -0.75 -0.744 ...
# $ selev2 : num [1:500] -0.793 -0.791 -0.782 -0.78 -0.774 ...

# Write JAGS model file
cat(file="model9.txt", "
model {
  # Priors and linear models: shared for models of both data sets
  alpha ~ dunif(-10, 10)                            # Abundance intercept
  mean.lam <- exp(alpha)
  beta ~ dnorm(0, 0.01)

  # Likelihoods for Data Sets 1 and 2
  # Note identical alpha and beta for both data sets
  for (i in 1:nsites1){                             # Data Set 1
  C[i] ~ dpois(lambda1[i])
    log(lambda1[i]) <- alpha + beta * selev1[i]
  }
  for (j in 1:nsites2){                             # Data Set 2
    y[j] ~ dbern(psi[j])
    cloglog(psi[j]) <- alpha + beta * selev2[j]

    # Alternative implementation of same model for Data Set 2
    # y[j] ~ dbern(psi[j])
    # psi[j] <- 1 - exp(-lambda2[j])
    # log(lambda2[j]) <- alpha + beta * selev2[j]
  }
}
")

# Initial values
inits <- function(){list(alpha=runif(1), beta=rnorm(1))}

# Parameters monitored
parameters <- c("mean.lam", "alpha", "beta")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 10; na <- 1000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out10 <- jags(jags.data, inits, parameters, "model9.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out10) # Not shown
print(out10, 3)
#               mean    sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# mean.lam     1.902 0.079    1.747    1.901    2.061    FALSE 1 1.000  3000
# alpha        0.642 0.041    0.558    0.643    0.723    FALSE 1 1.000  3000
# beta        -2.127 0.079   -2.281   -2.127   -1.971    FALSE 1 1.001  2429


# Definition of negative log-likelihood (NLL) for Poisson GLM
NLL1 <- function(param, y, x){
  alpha <- param[1]                                 # Organize parameters in a vector
  beta <- param[2]
  lambda <- exp(alpha + beta * x)                   # Linear predictor
  L <- dpois(y, lambda)                             # Likelihood contribution for each datum
  return(-sum(log(L)))                              # NLL for all observations in data set
}

# Minimize NLL
inits <- c(alpha=0, beta=0)                         # Need to provide initial values
sol1 <- optim(inits, NLL1, y=C1, x=selev1, hessian=TRUE)
mle1 <- sol1$par                                    # Grab the MLEs
VC1 <- solve(sol1$hessian)                          # Get variance-covariance matrix
ASE1 <- sqrt(diag(VC1))                             # Extract asymptotic SEs
print(cbind(mle1, ASE1), 3)                         # Print MLEs and asymptotic SEs
#       mle1   ASE1
# alpha 0.62 0.0621
# beta -2.19 0.1169

# ~~~~ Compare with 'official' solution ~~~~
summary(glm(C1 ~ selev1, family=poisson(link="log")))

# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.62035    0.06213   9.985   <2e-16 ***
# selev1      -2.19121    0.11691 -18.742   <2e-16 ***
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Definition of negative log-likelihood for Bernoulli GLM with cloglog link
NLL2 <- function(param, y, x){
  alpha <- param[1]
  beta <- param[2]
  lambda <- exp(alpha + beta * x)                   # Linear predictor
  psi <- 1-exp(-lambda)                             # cloglog link
  L <- dbinom(y, 1, psi)                            # Likelihood contribution for each datum
  return(-sum(log(L)))                              # NLL for all observations in data set
}

# Minimize NLL
inits <- c(alpha=0, beta=0)
sol2 <- optim(inits, NLL2, y=y, x=selev2, hessian=TRUE)
mle2 <- sol2$par                                    # Grab MLEs
VC2 <- solve(sol2$hessian)                          # Get variance-covariance matrix
ASE2 <- sqrt(diag(VC2))                             # Extract asymptotic SEs
print(cbind(mle2, ASE2), 3)                         # Print MLEs and asymptotic SEs
#        mle2   ASE2
# alpha 0.598 0.0794
# beta -1.953 0.1765

# ~~~~ Compare with 'official' solution ~~~~
summary(glm(y ~ selev2, family=binomial(link="cloglog")))

# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.59805    0.07968   7.506  6.1e-14 ***
# selev2      -1.95305    0.17628 -11.079  < 2e-16 ***
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Definition of the joint NLL for the integrated model
NLLjoint <- function(param, y1, x1, y2, x2){
  # Definition of elements in param vector (shared between data sets)
  alpha <- param[1]
  beta <- param[2]

  # Likelihood for the Poisson GLM for data set 1 (y1, x1)
  lambda1 <- exp(alpha + beta * x1)
  L1 <- dpois(y1, lambda1)                          # Likelihood contribution each datum

  # Likelihood for the cloglog Bernoulli GLM for data set 2 (y2, x2)
  lambda2 <- exp(alpha + beta * x2)
  psi <- 1-exp(-lambda2)
  L2 <- dbinom(y2, 1, psi)                          # Likelihood contribution each datum

  # Joint log-likelihood and joint NLL: here you can see that sum!
  JointLL <- sum(log(L1)) + sum(log(L2))            # Joint LL
  return(-JointLL)                                  # Return joint NLL
}

# Minimize NLLjoint
inits <- c(alpha=0, beta=0)
solJoint <- optim(inits, NLLjoint, y1=C1, x1=selev1, y2=y, x2=selev2, hessian=TRUE)

# Get MLE and asymptotic SE and print them
mleJoint <- solJoint$par                            # Grab MLEs
VC <- solve(solJoint$hessian)                       # Get variance-covariance matrix
ASE <- sqrt(diag(VC))                               # Extract asymptotic SEs
print(cbind(mleJoint, ASE), 3)                      # Print MLEs and asymptotic SEs
#       mleJoint    ASE
# alpha    0.642 0.0408
# beta    -2.130 0.0804

# ~~~~ Bayesian estimates for comparison ~~~~
#               mean    sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# mean.lam     1.902 0.079    1.747    1.901    2.061    FALSE 1 1.000  3000
# alpha.lam    0.642 0.041    0.558    0.643    0.723    FALSE 1 1.000  3000
# beta.lam    -2.127 0.079   -2.281   -2.127   -1.971    FALSE 1 1.001  2429
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
