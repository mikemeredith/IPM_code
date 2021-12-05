# Schaub & Kéry (2022) Integrated Population Models
# Chapter 8 : Integrated population models with density-dependence
# ----------------------------------------------------------------

# Run time for test script 12 mins, full run 1.7 hrs

library(IPMbook) ; library(jagsUI)

# ~~~ this uses data from section 8.2.2 ~~~
data(redbacked)
jags.data <- with(redbacked, list(n.occasions=ncol(marr.a), marr.j=marr.j,
    marr.a=marr.a, rel.j=rowSums(marr.j), rel.a=rowSums(marr.a),
    C=count, J=J, B=B, pNinit=dUnif(1, 50), mean.C=mean(count)))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 8.2 Density-dependence in red-backed shrikes
# ============================================

# 8.2.4 Modeling density-dependence in immigration
# ------------------------------------------------

# Write JAGS model file
cat(file="model2.txt", "
model {
  # Priors and linear models
  # Models for demographic rates
  for (t in 1:(n.occasions-1)){
    logit(phij[t]) <- alpha[1] + beta[1] * (N[t] - mean.C) + err.phij[t]
    err.phij[t] ~ dnorm(0, tau.phij)
    logit(phia[t]) <- alpha[2] + beta[2] * (N[t] - mean.C) + err.phia[t]
    err.phia[t] ~ dnorm(0, tau.phia)
    logit(pj[t]) <- logit.b0.pj + err.pj[t]
    err.pj[t] ~ dnorm(0, tau.pj)
    logit(pa[t]) <- logit.b0.pa + err.pa[t]
    err.pa[t] ~ dnorm(0, tau.pa)
  }
  for (t in 1:n.occasions){
    log(f[t]) <- alpha[3] + beta[3] * (N[t] - mean.C) + err.f[t]
    err.f[t] ~ dnorm(0, tau.f)
    log(omega[t]) <- alpha[4] + beta[4] * log(S[t] + R[t]) + err.om[t]
    err.om[t] ~ dnorm(0, tau.om)
  }

  # Priors for variance parameters (hyperparameters)
  tau.phij <- pow(sigma.phij, -2)
  sigma.phij ~ dunif(0, 10)
  sigma2.phij <- pow(sigma.phij, 2)
  tau.phia <- pow(sigma.phia, -2)
  sigma.phia ~ dunif(0, 10)
  sigma2.phia <- pow(sigma.phia, 2)
  tau.f <- pow(sigma.f, -2)
  sigma.f ~ dunif(0, 10)
  sigma2.f <- pow(sigma.f, 2)
  tau.om <- pow(sigma.om, -2)
  sigma.om ~ dunif(0, 10)
  sigma2.om <- pow(sigma.om, 2)
  tau.pj <- pow(sigma.pj, -2)
  sigma.pj ~ dunif(0, 10)
  sigma2.pj <- pow(sigma.pj, 2)
  tau.pa <- pow(sigma.pa, -2)
  sigma.pa ~ dunif(0, 10)
  sigma2.pa <- pow(sigma.pa, 2)
  tau ~ dgamma(0.001, 0.001)
  sigma2 <- 1 / tau

  # Priors for the mean of resighting rates
  b0.pj ~ dunif(0, 1)
  b0.pa ~ dunif(0, 1)

  # Priors for regression parameters of the demographic rates
  for (j in 1:4){
    alpha[j] ~ dnorm(0, 0.001)
    beta[j] ~ dunif(-5, 5)
  }

  # Back-transformations
  logit.b0.pj <- logit(b0.pj)
  logit.b0.pa <- logit(b0.pa)

  # Population count data (state-space model)
  # Model for initial stage-spec. population sizes: discrete uniform priors
  R[1] ~ dcat(pNinit)                             # Local recruits
  S[1] ~ dcat(pNinit)                             # Surviving adults
  I[1] ~ dpois(omega[1])                          # Immigrants

  # Process model over time: our model of population dynamics
  for (t in 2:n.occasions){
    R[t] ~ dpois(f[t-1]/2 * phij[t-1] * N[t-1])   # No. local recruits
    S[t] ~ dbin(phia[t-1], N[t-1])                # No. surviving adults
    I[t] ~ dpois(omega[t])                        # No. immigrants
  }

  # Observation process
  for (t in 1:n.occasions){
    N[t] <- S[t] + R[t] + I[t]                    # Total number of breeding females
    logN[t] <- log(N[t])
    C[t] ~ dlnorm(logN[t], tau)
  }

  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    qj[t] <- 1-pj[t]                              # Probability of non-recapture (juv)
    qa[t] <- 1-pa[t]                              # Probability of non-recapture (ad)
    pr.j[t,t] <- phij[t] * pj[t]
    pr.a[t,t] <- phia[t] * pa[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- phij[t] * prod(phia[(t+1):j]) * qj[t] * prod(qa[t:(j-1)]) * pa[j] / qa[t]
      pr.a[t,j] <- prod(phia[t:j]) * prod(qa[t:(j-1)]) * pa[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  }

  # Productivity data (Poisson regression model)
  for (t in 1:n.occasions){
    J[t] ~ dpois(B[t] * f[t])
  }
}
")

# Define initial values
inits <- function() {list(b0.pj=runif(1,0.2,0.7), b0.pa=runif(1,0.1,0.6), sigma.pa=runif(1,0,1),
    I=rep(10, 36))}

# Define parameters to be monitored
parameters <- c("phij", "phia", "f", "omega", "b0.pj", "b0.pa", "sigma2.phij", "sigma2.phia", "sigma2.f",
    "sigma2.om", "sigma2.pj", "sigma2.pa", "R", "S", "I", "N", "pa", "pj", "sigma2", "alpha", "beta")

# MCMC settings
# ni <- 110000; nb <- 10000; nc <- 3; nt <- 20; na <- 1000
ni <- 11000; nb <- 1000; nc <- 3; nt <- 2; na <- 1000  # ~~~ for testing, 11 mins

# Call JAGS (ART 139 min), check convergence and summarize posteriors
out2 <- jags(jags.data, inits, parameters, "model2.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out2)

print(out2, 3)

#                mean     sd    2.5%     50%    97.5% overlap0     f  Rhat n.eff
# phij[1]       0.065  0.025   0.026   0.061    0.126    FALSE 1.000 1.000 15000
# phij[2]       0.072  0.027   0.031   0.068    0.137    FALSE 1.000 1.000  7918
# ... [output truncated] ...
# phij[34]      0.062  0.021   0.028   0.060    0.110    FALSE 1.000 1.001  3062
# phij[35]      0.062  0.023   0.025   0.059    0.116    FALSE 1.000 1.000  6145
# phia[1]       0.369  0.057   0.259   0.368    0.488    FALSE 1.000 1.000  4193
# phia[2]       0.326  0.054   0.218   0.328    0.430    FALSE 1.000 1.000  8420
# ... [output truncated] ...
# phia[34]      0.395  0.052   0.298   0.393    0.506    FALSE 1.000 1.000 15000
# phia[35]      0.384  0.055   0.280   0.381    0.499    FALSE 1.000 1.001  2657
# fec[1]        3.696  0.330   3.083   3.687    4.370    FALSE 1.000 1.000 14422
# fec[2]        2.523  0.217   2.121   2.517    2.962    FALSE 1.000 1.000 15000
# ... [output truncated] ...
# fec[35]       3.321  0.231   2.886   3.314    3.787    FALSE 1.000 1.000 15000
# fec[36]       3.065  0.229   2.637   3.056    3.533    FALSE 1.000 1.000 15000
# omega[1]     27.684  2.988  21.668  27.717   33.504    FALSE 1.000 1.000 15000
# omega[2]     30.750  3.314  25.235  30.428   38.259    FALSE 1.000 1.000 10076
# ... [output truncated] ...
# omega[35]    32.034  3.666  25.190  31.864   39.659    FALSE 1.000 1.003   806
# omega[36]    30.940  3.413  24.122  30.908   37.913    FALSE 1.000 1.004   580
# b0.pj         0.442  0.069   0.316   0.438    0.588    FALSE 1.000 1.000 11575
# b0.pa         0.645  0.039   0.569   0.645    0.721    FALSE 1.000 1.000  6175
# sigma2.phij   0.240  0.179   0.004   0.206    0.684    FALSE 1.000 1.002  1189
# sigma2.phia   0.097  0.066   0.006   0.084    0.259    FALSE 1.000 1.001  3182
# sigma2.fec    0.051  0.015   0.028   0.048    0.087    FALSE 1.000 1.000 15000
# sigma2.om     0.009  0.013   0.000   0.004    0.044    FALSE 1.000 1.004  5289
# sigma2.pj     0.351  0.533   0.001   0.170    1.703    FALSE 1.000 1.000  6434
# sigma2.pa     0.077  0.107   0.000   0.039    0.384    FALSE 1.000 1.009   466
# R[1]         10.800  6.582   1.000  10.000   24.000    FALSE 1.000 1.000  6269
# R[2]          6.738  3.352   1.000   6.000   14.000    FALSE 1.000 1.000 13740
# ... [output truncated] ...
# R[35]         7.528  3.475   2.000   7.000   15.000    FALSE 1.000 1.001  2581
# R[36]         6.673  3.363   1.000   6.000   14.000    FALSE 1.000 1.000  7694
# S[1]         10.708  6.560   1.000  10.000   24.000    FALSE 1.000 1.000 15000
# S[2]         19.683  3.902  12.000  20.000   27.000    FALSE 1.000 1.000 10316
# ... [output truncated] ...
# S[35]        26.945  4.501  18.000  27.000   36.000    FALSE 1.000 1.000 10698
# S[36]        25.381  4.504  17.000  25.000   35.000    FALSE 1.000 1.000  6205
# I[1]         26.254  4.517  17.000  26.000   35.000    FALSE 1.000 1.001  5823
# I[2]         35.011  4.432  26.000  35.000   44.000    FALSE 1.000 1.000 15000
# ... [output truncated] ...
# I[35]        32.919  4.739  23.000  33.000   42.000    FALSE 1.000 1.001  2209
# I[36]        30.141  4.559  21.000  30.000   39.000    FALSE 1.000 1.001  2355
# N[1]         47.762  2.400  44.000  48.000   53.000    FALSE 1.000 1.000 15000
# N[2]         61.431  2.750  56.000  62.000   67.000    FALSE 1.000 1.000 15000
# ... [output truncated] ...
# N[35]        67.391  3.000  61.000  68.000   73.000    FALSE 1.000 1.000 11563
# N[36]        62.195  2.836  57.000  62.000   68.000    FALSE 1.000 1.000 15000
# pa[1]         0.633  0.070   0.475   0.637    0.765    FALSE 1.000 1.003   801
# pa[2]         0.633  0.072   0.470   0.638    0.766    FALSE 1.000 1.002  4871
# ... [output truncated] ...
# pa[34]        0.672  0.068   0.553   0.664    0.830    FALSE 1.000 1.006   899
# pa[35]        0.651  0.065   0.522   0.649    0.788    FALSE 1.000 1.001  3339
# pj[1]         0.469  0.133   0.238   0.454    0.793    FALSE 1.000 1.000 15000
# pj[2]         0.437  0.119   0.205   0.433    0.703    FALSE 1.000 1.001 15000
# ... [output truncated] ...
# pj[34]        0.441  0.117   0.217   0.436    0.696    FALSE 1.000 1.001  4079
# pj[35]        0.454  0.127   0.223   0.443    0.750    FALSE 1.000 1.000 11030
# sigma2.c      0.002  0.002   0.000   0.002    0.009    FALSE 1.000 1.007  6627
# alpha[1]     -2.833  0.181  -3.204  -2.827   -2.494    FALSE 1.000 1.001  2407
# alpha[2]     -0.540  0.091  -0.720  -0.540   -0.365    FALSE 1.000 1.000 12706
# alpha[3]      1.014  0.041   0.932   1.015    1.094    FALSE 1.000 1.000  8618
# alpha[4]      2.488  0.628   1.091   2.559    3.523    FALSE 1.000 1.021   137
# beta[1]       0.006  0.013  -0.020   0.006    0.032     TRUE 0.667 1.001  5265
# beta[2]       0.002  0.007  -0.012   0.002    0.017     TRUE 0.614 1.002   924
# beta[3]      -0.003  0.004  -0.010  -0.003    0.004     TRUE 0.799 1.000  8207
# beta[4]       0.275  0.197  -0.056   0.253    0.709     TRUE 0.939 1.021   136

# ~~~~ code for Figure 8.5 ~~~~
n <- length(out2$mean$N)
# Compute the number of locals
loc <- out2$sims.list$S + out2$sims.list$R
lower <- upper <- m <- numeric()
for (t in 1:n){
  lower[t] <- quantile(loc[,t], 0.025)
  upper[t] <- quantile(loc[,t], 0.975)
  m[t] <- mean(loc[,t])
}

library(scales)
library(plotrix)
op <- par(mfrow=c(1,2), mar=c(4.5,4,2,1), cex=1.15)
plot(density(out2$sims.list$beta[,4]), type="l", lwd=2, xlim=c(-0.5, 1), main=NA,
    xlab=expression('Strength of density-dependence ('*beta[I]*')'),
    ylab="Density", col="black", axes=FALSE)
axis(1)
axis(2, las=1)
corner.label("A", font=2, cex=1.15)

nx <- min(lower):max(upper)
I.pred <- matrix(NA, nrow=out2$mcmc.info$n.samples, ncol=length(nx))
for (i in 1:length(nx)){
  I.pred[,i] <- exp(out2$sims.list$alpha[,4] + out2$sims.list$beta[,4] * log(nx[i]))
}
quant <- function (x) quantile(x, c(0.025, 0.975))

plot(NA, ylab="Number of immigrants", xlab="Number of locals", axes=FALSE,
    ylim=c(10, max(out2$q97.5$I)), xlim=c(min(lower), max(upper)))
polygon(x=c(nx, max(upper):min(lower)), y=c(apply(I.pred, 2, quant)[1,],
    apply(I.pred, 2, quant)[2,length(nx):1]), border=NA, col=alpha("salmon2", 0.3))
lines(y=apply(I.pred, 2, mean), x=nx, lwd=2, col="salmon2")
segments(lower, out2$mean$I[1:n], upper, out2$mean$I[1:n], col="grey50")
segments(m, out2$q2.5$I[1:n], m, out2$q97.5$I[1:n], col="grey50")
points(y=out2$mean$I, x=m, pch=16, cex=1.3)
axis(2, las=1)
axis(1)
corner.label("B", font=2, cex=1.15)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
