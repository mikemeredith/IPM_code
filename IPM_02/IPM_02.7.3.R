# Schaub & Kéry (2022) Integrated Population Models
# Chapter 2 : Bayesian statistical modeling using JAGS
# ----------------------------------------------------

# Run time approx. 2 mins

library(IPMbook) ; library(jagsUI)

# 2.7 Using JAGS to fit simple statistical models from R: GLMs and GLMMs
# ======================================================================

# 2.7.3 Binomial generalized linear models
# ----------------------------------------

# Simulate a data set of binomial counts
set.seed(99)
nsites <- 1000                            # Number of sites

# Simulate true states
mean.abundance <- 6                                 # Mean abundance across sites
N <- rpois(nsites, mean.abundance)                  # Realized abundance N is Poisson
table(N)                                            # Note 7 unoccupied sites
# N
# 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  21
# 7  13  38 109 141 157 144 126 101  69  47  24  14   7   1   1   1

# Simulate one count per site between March and October
# Simulate survey dates
survey.date <- sort(round(runif(nsites, 60, 300)))  # Survey dates
sdate <- as.numeric(scale(survey.date))             # Scaled survey date (linear)
sdate2 <- as.numeric(sdate^2)                       # Scaled survey date (squared)

# Simulate seasonal profile of p
mean.p <- 0.05                                      # Detection probability in the middle of the season
beta1.lp <- -0.3                                    # Linear effect of scaled survey date
beta2.lp <- 0.9                                     # Quadratic effect of scaled survey date
p <- plogis(qlogis(mean.p) + beta1.lp * sdate + beta2.lp * sdate2)

# Simulate one binomial count at each site
C <- rbinom(nsites, N, p)

# ~~~~ extra code for Figure 2.15 ~~~~
# Plot abundance, detection, seasonal profile of counts
#  and counts versus abundance
op <- par(mfrow = c(2, 2), mar=c(5, 5, 2, 2), cex.lab=1.5, cex.axis=1.5,
    cex.main=1.5, las=1)
plot(survey.date, N, pch=16, col=rgb(0, 0, 0, 0.4), cex=1.5, xlab="Survey date",
    ylab=expression(paste('Snake abundance (', italic(N), ')')), frame=FALSE,
    ylim=c(0, max(N)))
plot(survey.date, p, type="l", lty=1, lwd=3, col="red", xlab="Survey date",
    ylab=expression(paste('Detection prob. (', italic(p), ')')), frame=FALSE,
    ylim=c(0, 0.6))
plot(survey.date, C, pch=16, col=rgb(0, 0, 0, 0.4), cex=1.5, xlab="Survey date",
    ylab=expression(paste('Snake count (', italic(C), ')')), frame=FALSE,
    ylim=c(0, max(C)))
plot(jitter(N), C, pch=16, col=rgb(0, 0, 0, 0.4), cex=1.5,
    xlab=expression(paste('Snake abundance (', italic(N), ')')),
    ylab=expression(paste('Snake count (', italic(C), ')')), frame=FALSE,
    xlim=range(N), ylim=range(N))
abline(0, 1, lwd=2, col="red")
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Bundle data
jags.data <- list(C=C, N=N, sdate=sdate, sdate2=sdate2, n=length(C))
str(jags.data)
# List of 5
# $ C     : int [1:1000] 3 0 3 9 3 2 4 3 3 1 ...
# $ N     : int [1:1000] 6 3 7 13 6 11 7 5 5 4 ...
# $ sdate : num [1:1000] -1.72 -1.72 -1.7 -1.7 -1.7 ...
# $ sdate2: num [1:1000] 2.96 2.96 2.91 2.91 2.91 ...
# $ n     : int 1000

# Write JAGS model file
cat(file="model3.txt", "
model {
  # Priors and linear models
  alpha <- logit(mean.p)
  mean.p ~ dunif(0, 1)                    # logit-linear intercept in prob. scale
  beta[1] ~ dunif(-10, 10)                # logit-linear effect of date (linear)
  beta[2] ~ dunif(-10, 10)                # logit-linear effect of date (squared)

  # Likelihood of the binomial GLM
  for (i in 1:n){
  C[i] ~ dbinom(p[i], N[i])
    logit(p[i]) <- alpha + beta[1] * sdate[i] + beta[2] * sdate2[i]
    # p[i] <- ilogit(alpha + beta[1] * sdate[i] + beta[2] * sdate2[i]) # same
  }
}
")

# Initial values
inits <- function(){list(mean.p=runif(1, 0, 0.3), beta=runif(2, 0, 2))}

# Parameters monitored
parameters <- c("alpha", "beta", "mean.p", "p")

# MCMC settings
ni <- 20000; nb <- 10000; nc <- 3; nt <- 2; na <- 1000

# Call JAGS from R (ART 1 min), check convergence and summarize posteriors
out4 <- jags(jags.data, inits, parameters, "model3.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out4) # Not shown (and may have to do jagsUI::traceplot)
print(out4, 3)

#             mean    sd     2.5%      50%    97.5% overlap0 f Rhat n.eff
# alpha      -2.978 0.073   -3.123   -2.977   -2.833    FALSE 1    1 13943
# beta[1]    -0.304 0.032   -0.367   -0.304   -0.240    FALSE 1    1 15000
# beta[2]     0.950 0.042    0.868    0.950    1.032    FALSE 1    1 11003
# mean.p      0.049 0.003    0.042    0.048    0.056    FALSE 1    1 14542
# p[1]        0.587 0.021    0.544    0.588    0.629    FALSE 1    1  9071
# p[2]        0.587 0.021    0.544    0.588    0.629    FALSE 1    1  9071
# p[3]        0.575 0.021    0.532    0.575    0.616    FALSE 1    1  9133
# [... output truncated ...]


# Definition of NLL for a logistic regression with two covariates
NLL <- function(param, y, N, x1, x2){
  alpha <- param[1]
  beta1 <- param[2]
  beta2 <- param[3]
  psi <- plogis(alpha + beta1 * x1 + beta2 * x2)
  L <- dbinom(y, N, psi)                    # Likelihood contribution for 1 datum
  NLL <- -sum(log(L))                       # NLL for all data points (here, 1000)
  return(NLL)
}

# Minimize NLL to obtain the MLEs and their asymptotic SEs
inits <- c(alpha=0, beta1=0, beta2=0)
sol <- optim(inits, NLL, y=C, N=N, x1=sdate, x2=sdate2, hessian=TRUE)
mle <- sol$par                              # Grab MLEs
VC <- solve(sol$hessian)                    # Get variance-covariance matrix
ASE <- sqrt(diag(VC))                       # Extract asymptotic SEs
print(cbind(mle, ASE), 3)                   # Print MLEs and SEs

#          mle    ASE
# alpha -2.979 0.0745
# beta1 -0.304 0.0323
# beta2  0.951 0.0423

# Get frequentist MLEs for comparison
summary(glm(cbind(C, N-C) ~ sdate + sdate2, family='binomial'))

# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept) -2.97874    0.07448 -39.991   <2e-16 ***
# sdate       -0.30354    0.03227  -9.407   <2e-16 ***
# sdate2       0.95098    0.04229  22.487   <2e-16 ***

# Plot seasonal detection profile: truth (red) and estimate (blue) (not shown)
plot(survey.date, p, type="l", lty=1, lwd=5, col="red", xlab="Survey date (day of year)",
    ylab="Detection probabilty", frame=FALSE, ylim=c(0, 0.6))
polygon(c(survey.date, rev(survey.date)), c(out4$q2.5$p, rev(out4$q97.5$p)), col="grey", border=NA)
lines(survey.date, p, type="l", lty=1, lwd=5, col="red")
lines(survey.date, out4$mean$p, type="l", lty=2, lwd=5, col="blue")
