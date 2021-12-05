# Schaub & Kéry (2022) Integrated Population Models
# Chapter 9 : Retrospective population analyses
# ---------------------------------------------

# Run time approx. 1 min

# 9.2 Correlations between demographic rates and population growth
# ================================================================

# Load the hoopoe data and produce data overview
library(IPMbook); library(jagsUI)
data(hoopoe)
str(hoopoe)
# List of 5
# $ ch       : int [1:3844, 1:16] 0 0 0 0 0 0 0 0 0 0 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:16] "2002" "2003" "2004" "2005" ...
# $ age      : int [1:3844] 1 1 1 1 1 1 1 1 1 1 ...
# $ count    : num [1:16] 34 46 68 93 88 87 85 78 82 84 ...
# $ reproAgg:List of 4
# ..$ J1: num [1:16] 73 188 243 320 261 222 206 154 278 220 ...
# ..$ J2: num [1:16] 51 83 101 182 256 226 206 278 237 244 ...
# ..$ B1: num [1:16] 15 23 36 50 44 37 39 30 48 48 ...
# ..$ B2: num [1:16] 6 13 12 23 38 38 32 41 39 41 ...
# $ reproInd:List of 3
# ..$ f : int [1:1092] 6 11 11 3 15 6 11 12 7 6 ...
# ..$ id : num [1:1092] 483 486 530 531 532 533 534 535 536 538 ...
# ..$ year: int [1:1092] 2002 2002 2002 2002 2002 2002 2002 2002 2002 ...

# Produce age-dependent m-arrays
marr <- marrayAge(hoopoe$ch, hoopoe$age, 2)
marr.j <- marr[,,1]
marr.a <- marr[,,2]

# Bundle data
jags.data <- with(hoopoe$reproAgg, list(marr.j=marr.j, marr.a=marr.a, n.occasions=ncol(marr.j),
    rel.j=rowSums(marr.j), rel.a=rowSums(marr.a), J1=J1, B1=B1, J2=J2, B2=B2, count=hoopoe$count,
    pNinit=dUnif(1, 50)))
str(jags.data)                                    # Remind ourselves of how the data look like

# Write JAGS model file
cat(file="model1.txt", "
model {
  # Priors and linear models
  # For mean and variance parameters
  mean.logit.phij <- logit(mean.phij)
  mean.phij ~ dunif(0, 1)
  mean.logit.phia <- logit(mean.phia)
  mean.phia ~ dunif(0, 1)
  mean.log.f1 <- log(mean.f1)
  mean.f1 ~ dunif(0, 10)
  mean.log.f2 <- log(mean.f2)
  mean.f2 ~ dunif(0, 10)
  mean.log.omega <- log(mean.omega)
  mean.omega ~ dunif(0, 3)
  mean.p ~ dunif(0, 1)
  mean.pj ~ dunif(0, 1)
  sigma.phij ~ dunif(0, 3)
  tau.phij <- pow(sigma.phij, -2)
  sigma.phia ~ dunif(0, 3)
  tau.phia <- pow(sigma.phia, -2)
  sigma.omega ~ dunif(0, 3)
  tau.omega <- pow(sigma.omega, -2)
  sigma.f1 ~ dunif(0, 3)
  tau.f1 <- pow(sigma.f1, -2)
  sigma.f2 ~ dunif(0, 3)
  tau.f2 <- pow(sigma.f2, -2)
  sigma ~ dunif(0.4, 20)                          # lower bound >0
  tau <- pow(sigma, -2)

  # Linear models for demographic parameters
  for (t in 1:(n.occasions-1)){
    logit.phij[t] ~ dnorm(mean.logit.phij, tau.phij)
    phij[t] <- ilogit(logit.phij[t])
    logit.phia[t] ~ dnorm(mean.logit.phia, tau.phia)
    phia[t] <- ilogit(logit.phia[t])
    log.omega[t] ~ dnorm(mean.log.omega, tau.omega)
    omega[t] <- exp(log.omega[t])
    p[t] <- mean.p
    pj[t] <- mean.pj
  }
  for (t in 1:n.occasions){
    log.f1[t] ~ dnorm(mean.log.f1, tau.f1)
    f1[t] <- exp(log.f1[t])
    log.f2[t] ~ dnorm(mean.log.f2, tau.f2)
    f2[t] <- exp(log.f2[t])
  }

  # Population count data (state-space model)
  # Model for initial stage-spec. population sizes: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)
  N[3,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois(f1[t] / 2 * phij[t] * (N[1,t] + N[3,t]) + f2[t] / 2 * phij[t] * N[2,t])
    N[2,t+1] ~ dbin(phia[t], (N[1,t] + N[2,t] + N[3,t]))
    N[3,t+1] ~ dpois((N[1,t] + N[2,t] + N[3,t]) * omega[t])
  }

  # Observation model
  for (t in 1:n.occasions){
    count[t] ~ dnorm(Ntot[t], tau)
    Ntot[t] <- N[1,t] + N[2,t] + N[3,t]           # total population size
  }

  # Productivity data (Poisson regression model)
  for (t in 1:n.occasions){
    J1[t] ~ dpois(f1[t] * B1[t])
    J2[t] ~ dpois(f2[t] * B2[t])
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
    qj[t] <- 1 - pj[t]
    q[t] <- 1 - p[t]
    pr.j[t,t] <- phij[t] * pj[t]
    pr.a[t,t] <- phia[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- phij[t] * prod(phia[(t+1):j]) * qj[t] * prod(q[t:(j-1)]) * p[j] / q[t]
      pr.a[t,j] <- prod(phia[t:j]) * prod(q[t:(j-1)]) * p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1 - sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1 - sum(pr.a[t,1:(n.occasions-1)])
  } #t

  # Derived parameters
  # Annual population growth rate
  for (t in 1:(n.occasions-1)){
    lambda[t] <- Ntot[t+1] / Ntot[t]
  }
}
")

# Initial values
inits <- function(){list(mean.phij=runif(1, 0.05, 0.2), mean.phia=runif(1, 0.3, 0.5))}

# Parameters monitored
parameters <- c("mean.phij", "sigma.phij", "mean.phia", "sigma.phia", "mean.omega", "sigma.omega",
    "mean.f1", "sigma.f1", "mean.f2", "sigma.f2", "mean.p", "mean.pj", "sigma", "phij", "phia", "omega",
    "f1", "f2", "N", "Ntot", "lambda")

# MCMC settings
ni <- 20000; nb <- 10000; nc <- 3; nt <- 4; na <- 2000

# Call JAGS (ART 1 min), check convergence and summarize posteriors
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1)
print(out1, 3)
#                mean     sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# mean.phij     0.123  0.012   0.100   0.123   0.148    FALSE 1 1.009   249
# sigma.phij    0.280  0.103   0.102   0.272   0.509    FALSE 1 1.001  3609
# mean.phia     0.364  0.020   0.325   0.364   0.406    FALSE 1 1.000  7500
# sigma.phia    0.175  0.107   0.012   0.163   0.419    FALSE 1 1.005   553
# mean.omega    0.281  0.047   0.187   0.282   0.372    FALSE 1 1.009   244
# sigma.omega   0.236  0.190   0.005   0.196   0.701    FALSE 1 1.009   440
# mean.f1       5.463  0.286   4.919   5.457   6.041    FALSE 1 1.000  5840
# sigma.f1      0.192  0.046   0.122   0.186   0.303    FALSE 1 1.002  4224
# mean.f2       6.137  0.381   5.422   6.125   6.946    FALSE 1 1.000  7500
# sigma.f2      0.226  0.054   0.144   0.218   0.354    FALSE 1 1.000  7500
# mean.p        0.794  0.033   0.727   0.795   0.856    FALSE 1 1.000  7500
# mean.pj       0.604  0.040   0.527   0.604   0.682    FALSE 1 1.003   629
# sigma         3.871  2.636   0.477   3.375  10.200    FALSE 1 1.001  2590
# phij[1]       0.124  0.025   0.079   0.122   0.177    FALSE 1 1.004   539
# [ ... output truncated ... ]
# phij[15]      0.108  0.023   0.066   0.108   0.154    FALSE 1 1.001  2463
# phia[1]       0.362  0.042   0.277   0.361   0.449    FALSE 1 1.001  7500
# [ ... output truncated ... ]
# phia[15]      0.403  0.050   0.331   0.394   0.524    FALSE 1 1.001  1945
# omega[1]      0.336  0.119   0.173   0.310   0.656    FALSE 1 1.001  1814
# [ ... output truncated ... ]
# omega[15]     0.284  0.077   0.138   0.281   0.458    FALSE 1 1.003  1776
# f1[1]         5.095  0.511   4.128   5.079   6.132    FALSE 1 1.000  7500
# [ ... output truncated ... ]
# f1[16]        5.206  0.415   4.440   5.187   6.054    FALSE 1 1.000  7500
# f2[1]         7.822  1.000   6.096   7.758   9.987    FALSE 1 1.000  4283
# [ ... output truncated ... ]
# f2[16]        6.231  0.434   5.404   6.215   7.116    FALSE 1 1.000  7500
# N[1,1]       10.957  8.469   1.000   9.000  31.000    FALSE 1 1.001  1481
# N[2,1]       15.212  9.165   1.000  14.000  33.000    FALSE 1 1.002  1205
# N[3,1]       11.794  8.591   1.000  10.000  31.000    FALSE 1 1.001  2304
# [ ... output truncated ... ]
# N[1,16]      14.032  4.215   6.000  14.000  23.000    FALSE 1 1.000  3704
# N[2,16]      20.046  3.843  13.000  20.000  28.000    FALSE 1 1.000  7500
# N[3,16]      14.007  4.527   6.000  14.000  23.525    FALSE 1 1.001  2786
# Ntot[1]      37.963  5.226  32.000  36.000  52.000    FALSE 1 1.000  3053
# [ ... output truncated ... ]
# Ntot[16]     48.085  3.864  40.000  48.000  57.000    FALSE 1 1.000  7500
# lambda[1]     1.252  0.143   0.937   1.278   1.500    FALSE 1 1.000  5035
# [ ... output truncated ... ]
# lambda[15]    0.976  0.087   0.793   0.980   1.150    FALSE 1 1.001  7500

# ~~~~ save output for use later ~~~~~
save(out1, file="ResultsHoopoe.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n.draws <- out1$mcmc.info$n.samples               # Determine number of MCMC draws
corr <- matrix(NA, ncol=5, nrow=n.draws)          # Create object to hold results
draws <- out1$sims.list                           # Dig out MCMC samples
for (s in 1:n.draws){ # Loop over all MCMC draws and get correlations
  corr[s,1] <- cor(draws$phij[s,], draws$lambda[s,])
  corr[s,2] <- cor(draws$phia[s,], draws$lambda[s,])
  corr[s,3] <- cor(draws$omega[s,], draws$lambda[s,])
  corr[s,4] <- cor(draws$f1[s,1:15], draws$lambda[s,])
  corr[s,5] <- cor(draws$f2[s,1:15], draws$lambda[s,])
}

# Calculate posterior summaries for the correlation coefficients
# Posterior means
apply(corr, 2, mean)
# 95% credible intervals
cri <- function(x) quantile(x, c(0.025, 0.975))
apply(corr, 2, cri)

# ~~~~ Calculation of correlation coefficients without uncertainty  ~~~~
cor(out1$mean$phij, out1$mean$lambda)
cor(out1$mean$phia, out1$mean$lambda)
cor(out1$mean$omega, out1$mean$lambda)
cor(out1$mean$f1[1:15], out1$mean$lambda)
cor(out1$mean$f2[1:15], out1$mean$lambda)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~ code for Figure 9.2 ~~~~
cri <- function(x) quantile(x, c(0.025, 0.975))
cri.lambda <- apply(out1$sims.list$lambda, 2, cri)

op <- par(mfrow=c(2, 3))
layout(matrix(1:6, 2, 3, byrow=TRUE), widths=c(1.05, 1, 1), heights=c(1, 1), TRUE)

cri.rate <- apply(out1$sims.list$phij, 2, cri)
plot(NA, ylim=range(cri.lambda), xlim=range(cri.rate),
    ylab=expression('Population growth rate ('*lambda*')'),
    xlab=expression('Juvenile survival ('*phi[italic(j)]*')'), axes=FALSE)
axis(1)
axis(1, at=c(0.075, 0.125, 0.175, 0.225), labels=NA, tcl=-0.25)
axis(2, las=1)
segments(out1$mean$phij, cri.lambda[1,], out1$mean$phij, cri.lambda[2,], col="grey60")
segments(cri.rate[1,], out1$mean$lambda, cri.rate[2,], out1$mean$lambda, col="grey60")
points(y=out1$mean$lambda, x=out1$mean$phij, pch=16)

cri.rate <- apply(out1$sims.list$phia, 2, cri)
plot(NA, ylim=range(cri.lambda), xlim=range(cri.rate), ylab=NA,
    xlab=expression('Adult survival ('*phi[italic(a)]*')'), axes=FALSE)
axis(1)
axis(2, las=1)
segments(out1$mean$phia, cri.lambda[1,], out1$mean$phia, cri.lambda[2,], col="grey60")
segments(cri.rate[1,], out1$mean$lambda, cri.rate[2,], out1$mean$lambda, col="grey60")
points(y=out1$mean$lambda, x=out1$mean$phia, pch=16)

cri.rate <- apply(out1$sims.list$omega, 2, cri)
plot(NA, ylim=range(cri.lambda), xlim=range(cri.rate),
    ylab=NA, xlab=expression('Immigration rate ('*omega*')'), axes=FALSE)
axis(1)
axis(2, las=1)
segments(out1$mean$omega, cri.lambda[1,], out1$mean$omega, cri.lambda[2,], col="grey60")
segments(cri.rate[1,], out1$mean$lambda, cri.rate[2,], out1$mean$lambda, col="grey60")
points(y=out1$mean$lambda, x=out1$mean$omega, pch=16)

cri.rate <- apply(out1$sims.list$f1, 2, cri)
plot(NA, ylim=range(cri.lambda), xlim=range(cri.rate),
    ylab=expression('Population growth rate ('*lambda*')'),
    xlab=expression('Productivity 1y ('*italic(f)[1]*')'), axes=FALSE)
axis(1)
axis(2, las=1)
segments(out1$mean$f1[1:15], cri.lambda[1,], out1$mean$f1[1:15], cri.lambda[2,], col="grey60")
segments(cri.rate[1,1:15], out1$mean$lambda, cri.rate[2,1:15], out1$mean$lambda, col="grey60")
points(y=out1$mean$lambda, x=out1$mean$f1[1:15], pch=16)

cri.rate <- apply(out1$sims.list$f2, 2, cri)
plot(NA, ylim=range(cri.lambda), xlim=range(cri.rate), ylab=NA,
    xlab=expression('Productivity ad ('*italic(f)[2]*')'), axes=FALSE)
axis(1)
axis(2, las=1)
segments(out1$mean$f2[1:15], cri.lambda[1,], out1$mean$f2[1:15], cri.lambda[2,], col="grey60")
segments(cri.rate[1,1:15], out1$mean$lambda, cri.rate[2,1:15], out1$mean$lambda, col="grey60")
points(y=out1$mean$lambda, x=out1$mean$f2[1:15], pch=16)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
