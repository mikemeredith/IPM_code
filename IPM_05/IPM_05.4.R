# Schaub & Kéry (2022) Integrated Population Models
# Chapter 5 : Introduction to integrated population models
# --------------------------------------------------------

# Run time approx. 1 min

library(IPMbook) ; library(jagsUI)

# ~~~ requires data prepared in section 5.2 ~~~
data(woodchat5)
marr <- marrayAge(woodchat5$ch, woodchat5$age)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5.4 The 3-step approach to integrated population modeling
# =========================================================

# 5.4.1 Development of a model that links demographic data with population size (no code)
# 5.4.2 Formulation of the likelihood for each available data set separately (no code)
# 5.4.3 Formulation of the joint likelihood (no code)

# 5.4.4 Writing the BUGS code for the Integrated Population Model
# ---------------------------------------------------------------

# Bundle data and produce data overview
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=woodchat5$repro[,1],
    year=woodchat5$repro[,2], age=woodchat5$repro[,3], C=woodchat5$count, pNinit=dUnif(1, 300))
str(jags.data)
# List of 10
# $ marr.j     : num [1:19, 1:20] 8 0 0 0 0 0 0 0 0 0 ...
# $ marr.a     : num [1:19, 1:20] 16 0 0 0 0 0 0 0 0 0 ...
# $ n.occasions: int 20
# $ rel.j      : num [1:19] 51 53 55 65 73 66 61 76 65 75 ...
# $ rel.a      : num [1:19] 36 39 44 61 61 50 43 61 51 53 ...
# $ J          : num [1:929] 6 2 2 5 3 5 3 2 3 2 ...
# $ year       : num [1:929] 1 1 1 1 1 1 1 1 1 1 ...
# $ age        : num [1:929] 1 1 1 1 1 1 1 1 1 1 ...
# $ C          : num [1:20] 91 119 131 88 139 145 148 116 112 106 ...
# $ pNinit     : num [1:300] 0.00333 0.00333 0.00333 0.00333 0.00333 ...

# Write JAGS model file
cat(file="model4.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.occasions-1)){
    logit.sj[t] ~ dnorm(mu.sj, tau.sj)
    sj[t] <- ilogit(logit.sj[t])          # Back-transformation from logit scale
    logit.sa[t] ~ dnorm(mu.sa, tau.sa)
    sa[t] <- ilogit(logit.sa[t])          # Back-transformation from logit scale
    p[t] <- mean.p
  }

  for (t in 1:n.occasions){
    log.f[1,t] ~ dnorm(mu.f[1], tau.f[1])
    f[1,t] <- exp(log.f[1,t])             # Back-transformation from log scale
    log.f[2,t] ~ dnorm(mu.f[2], tau.f[2])
    f[2,t] <- exp(log.f[2,t])             # Back-transformation from log scale
  }

  mean.sj ~ dunif(0, 1)
  mu.sj <- logit(mean.sj)                 # Logit transformation
  mean.sa ~ dunif(0, 1)
  mu.sa <- logit(mean.sa)                 # Logit transformation
  sigma.sj ~ dunif(0, 3)
  tau.sj <- pow(sigma.sj, -2)
  sigma.sa ~ dunif(0, 3)
  tau.sa <- pow(sigma.sa, -2)

  for (j in 1:2){
    mean.f[j] ~ dunif(0, 10)
    mu.f[j] <- log(mean.f[j])             # Log transformation
    sigma.f[j] ~ dunif(0, 3)
    tau.f[j] <- pow(sigma.f[j], -2)
  }

  mean.p ~ dunif(0, 1)

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for initial stage-spec. population sizes: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois(N[1,t] * f[1,t] / 2 * sj[t] + N[2,t] * f[2,t] / 2 * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Productivity data (Poisson regression model)
  for (i in 1:length(J)){
    J[i] ~ dpois(f[age[i],year[i]])
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
    q[t] <- 1 - p[t]                      # Probability of non-recapture
    pr.j[t,t] <- sj[t] * p[t]
    pr.a[t,t] <- sa[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t] * prod(sa[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(sa[t:j]) * prod(q[t:(j-1)]) * p[j]
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

  # Derived parameters
  # Annual population growth rate (added 0.001 to avoid possible division by 0)
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t] + 0.001)
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.f", "mean.p", "sigma.sj", "sigma.sa", "sigma.f",
    "sigma", "sj", "sa", "f", "N", "ann.growth.rate", "Ntot")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 2; na <- 1000

# Call JAGS (ART 1 min), check convergence and summarize posteriors
out4 <- jags(jags.data, inits, parameters, "model4.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out4)
print(out4, 3)
#                         mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# mean.sj                0.301  0.014    0.275    0.301    0.327    FALSE 1 1.001  2027
# mean.sa                0.542  0.016    0.511    0.542    0.573    FALSE 1 1.000  5642
# mean.f[1]              2.673  0.083    2.502    2.673    2.829    FALSE 1 1.006  1000
# mean.f[2]              3.675  0.110    3.468    3.672    3.900    FALSE 1 1.012   185
# mean.p                 0.605  0.019    0.567    0.605    0.643    FALSE 1 1.000 15000
# sigma.sj               0.114  0.077    0.007    0.103    0.291    FALSE 1 1.011   203
# sigma.sa               0.122  0.086    0.005    0.109    0.318    FALSE 1 1.001  4371
# sigma.f[1]             0.041  0.033    0.000    0.034    0.120    FALSE 1 1.001  2706
# sigma.f[2]             0.064  0.041    0.003    0.061    0.153    FALSE 1 1.002  1090
# sigma                 13.028  3.575    7.287   12.584   21.169    FALSE 1 1.001  1491
# sj[1]                  0.305  0.027    0.253    0.304    0.365    FALSE 1 1.005   980
# [ ... output truncated ... ]
# sj[19]                 0.304  0.025    0.256    0.303    0.360    FALSE 1 1.001 13466
# sa[1]                  0.554  0.036    0.490    0.548    0.641    FALSE 1 1.001  9699
# [ ... output truncated ... ]
# sa[19]                 0.538  0.032    0.468    0.539    0.604    FALSE 1 1.001  2810
# f[1,1]                 2.723  0.159    2.462    2.704    3.114    FALSE 1 1.002  5080
# f[2,1]                 4.012  0.400    3.513    3.905    5.014    FALSE 1 1.003  1260
# [ ... output truncated ... ]
# f[1,20]                2.691  0.140    2.435    2.683    3.010    FALSE 1 1.002  3466
# f[2,20]                3.814  0.229    3.442    3.783    4.342    FALSE 1 1.004   576
# N[1,1]                41.533 28.882    2.000   37.000  103.000    FALSE 1 1.002  1220
# N[2,1]                60.414 26.969    6.000   64.000  103.000    FALSE 1 1.002  1080
# [ ... output truncated ... ]
# N[1,20]               71.405  8.953   54.000   71.000   90.000    FALSE 1 1.000 15000
# N[2,20]               79.945  7.909   65.000   80.000   95.000    FALSE 1 1.001  5014
# ann.growth.rate[1]     1.130  0.106    0.936    1.124    1.351    FALSE 1 1.001  3874
# [ ... output truncated ... ]
# ann.growth.rate[19]    1.019  0.068    0.890    1.019    1.158    FALSE 1 1.000 15000
# Ntot[1]              101.947  9.925   84.000  101.000  123.000    FALSE 1 1.000 12582
# [ ... output truncated ... ]
# Ntot[20]             151.349 10.377  131.000  151.000  173.000    FALSE 1 1.000 15000


# ~~~~ code for Figure 5.4 ~~~~
mag <- 1.25
cex.tif <- mag * 1.25
lwd.tif <- 3 * mag
op <- par(mar=c(4, 4, 3, 0), las=1, cex=cex.tif, lwd=lwd.tif)
u <- col2rgb("grey82")
T <- length(woodchat5$count)
col.pol <- rgb(u[1], u[2], u[3], alpha=100, maxColorValue=255)
plot(out4$mean$Ntot, type="n",
    ylim=range(c(out4$q2.5$Ntot, out4$q97.5$Ntot, woodchat5$count)),
    ylab="Population size", xlab="Year", las=1, cex=1.5, axes=FALSE)
axis(2, las=1, lwd=lwd.tif)
axis(2, at=c(90, 110, 130, 150), labels=NA, tcl=-0.25, lwd=lwd.tif)
axis(1, at=1:T, labels=NA, tcl=-0.25, lwd=lwd.tif)
axis(1, at=c(5, 10, 15, 20), labels=c(5, 10, 15, 20), tcl=-0.5, lwd=lwd.tif)
polygon(c(1:T, T:1), c(out4$q2.5$Ntot, out4$q97.5$Ntot[T:1]), border=NA, col=col.pol)
points(out4$mean$Ntot, type="b", col="black", pch=16, lty=1, lwd=lwd.tif)
points(woodchat5$count, type="b", col="blue", pch=1, lty=2, lwd=lwd.tif)
legend("topleft", legend=c("Observed population counts", "Estimated population size"),
    pch=c(1, 16), lwd=c(lwd.tif, lwd.tif), col=c("blue", "black"), lty=c(2, 1), bty="n")
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
