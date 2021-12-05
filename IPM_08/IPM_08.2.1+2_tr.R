# Schaub & Kéry (2022) Integrated Population Models
# Chapter 8 : Integrated population models with density-dependence
# ----------------------------------------------------------------

# Run time for test script 4 mins, full run 30 mina

# 8.2 Density-dependence in red-backed shrikes
# ============================================

# 8.2.1 General population model (no code)

# 8.2.2 Modeling density-dependence in survival and productivity
# --------------------------------------------------------------

# Load the red-backed shrike data and produce data overview
library(IPMbook); library(jagsUI)
data(redbacked)
str(redbacked) # Not shown

# Bundle data
jags.data <- with(redbacked, list(n.occasions=ncol(marr.a), marr.j=marr.j, marr.a=marr.a,
    rel.j=rowSums(marr.j), rel.a=rowSums(marr.a), C=count, J=J, B=B, pNinit=dUnif(1, 50),
    mean.C=mean(count)))
str(jags.data)
# List of 10
# $ n.occasions: int 36
# $ marr.j     : num [1:35, 1:36] 2 0 0 0 0 0 0 0 0 0 ...
# $ marr.a     : num [1:35, 1:36] 6 0 0 0 0 0 0 0 0 0 ...
# $ rel.j      : num [1:35] 61 66 45 85 56 120 90 69 54 45 ...
# $ rel.a      : num [1:35] 35 39 17 34 32 47 45 46 37 27 ...
# $ C          : num [1:36] 47 62 62 64 72 74 68 65 38 45 ...
# $ J          : num [1:36] 112 114 88 167 108 223 171 127 102 84 ...
# $ B          : num [1:36] 29 46 50 56 65 62 63 60 33 41 ...
# $ pNinit     : num [1:50] 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 ...
# $ mean.C     : num 53.4

# Write JAGS model file
cat(file="model1.txt", "
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
    log(omega[t]) <- log.b0.om + err.om[t]
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

  # Priors for the mean of immigration and resighting rates
  b0.om ~ dunif(0, 75)
  b0.pj ~ dunif(0, 1)
  b0.pa ~ dunif(0, 1)

  # Priors for regression parameters of the demographic rates
  for (j in 1:3){
    alpha[j] ~ dnorm(0, 0.001)
    beta[j] ~ dunif(-5, 5)
  }

  # Back-transformations
  log.b0.om <- log(b0.om)
  logit.b0.pj <- logit(b0.pj)
  logit.b0.pa <- logit(b0.pa)

  # Population count data (state-space model)
  # Model for initial stage-spec. population sizes: discrete uniform priors
  R[1] ~ dcat(pNinit)                             # Local recruits
  S[1] ~ dcat(pNinit)                             # Surviving adults
  I[1] ~ dpois(omega[1])                          # Immigrants

  # Process model over time: our model of population dynamics
  for (t in 2:n.occasions){
    R[t] ~ dpois(f[t-1] / 2 * phij[t-1] * N[t-1]) # No. local recruits
    S[t] ~ dbin(phia[t-1], N[t-1])                # No. surviving adults
    I[t] ~ dpois(omega[t])                        # No. immigrants
  }

  # Observation model
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
      pr.j[t,j] <- phij[t] * prod(phia[(t+1):j]) * qj[t]*prod(qa[t:(j-1)]) * pa[j] / qa[t]
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
inits <- function() {list(b0.om=runif(1,0,1), b0.pj=runif(1,0.2,0.7), b0.pa=runif(1,0.1,0.6),
    sigma.pa=runif(1,0,1), I=rep(10, 36))}

# Define parameters to be monitored
parameters <- c("phij", "phia", "f", "omega", "b0.om", "b0.pj", "b0.pa", "sigma2.phij",
    "sigma2.phia", "sigma2.f", "sigma2.om", "sigma2.pj", "sigma2.pa", "R", "S", "I", "N", "pa", "pj",
    "sigma2", "alpha", "beta")

# MCMC settings
# ni <- 30000; nb <- 10000; nc <- 3; nt <- 4; na <- 1000
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000  # ~~~ for testing, 3.5 mins

# Call JAGS (ART 31 min), check convergence and produce figures
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1)

# Compute the probability that the betas are negative
mean(out1$sims.list$beta[,1] < 0)
# [1] 0.2459
mean(out1$sims.list$beta[,2] < 0)
# [1] 0.2258
mean(out1$sims.list$beta[,3] < 0)
# [1] 0.7961

# Calculate population growth rate
lam <- out1$sims.list$N[,-1] / out1$sims.list$N[,-ncol(out1$sims.list$N)]

# ~~~~ save output for use later ~~~~
save(out1, file="Shrike-mod1.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for figure 8.2 ~~~~
# load("Shrike-mod1.Rdata")  # if necessary

plot(density(out1$sims.list$beta[,1]), type="l", lwd=2, ylim=c(0, 120),
    xlim=c(-0.05, 0.05), main=NA,
    xlab=expression('Strength of density-dependence ('*beta*')'),
    ylab="Density", col="forestgreen", axes=FALSE)
lines(density(out1$sims.list$beta[,2]), type="l", lwd=2, col="red")
lines(density(out1$sims.list$beta[,3]), type="l", lwd=2, col="blue")
legend("topleft", lwd=rep(2,3), col=c("forestgreen", "red", "blue"),
    legend=c(expression('Juvenile survival ('*beta[1]*')'),
    expression('Adult survival ('*beta[2]*')'),
    expression('Productivity ('*beta[3]*')')), bty="n")
abline(v=0, lty=2)
axis(1)
axis(2, las=1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~ Predicted demographic rates ~~~~
n <- length(out1$mean$N)-1
nx.min <- round(min(out1$q2.5$N[1:n]))
nx.max <- round(max(out1$q97.5$N[1:n]))
nx <- nx.min:nx.max
phij.pred <- phia.pred <- fec.pred <- matrix(NA, nrow=out1$mcmc.info$n.samples, ncol=length(nx))
for (i in 1:length(nx)){
  phij.pred[,i] <- plogis(out1$sims.list$alpha[,1] + out1$sims.list$beta[,1] *
      (nx[i] -jags.data$mean.C))
  phia.pred[,i] <- plogis(out1$sims.list$alpha[,2] + out1$sims.list$beta[,2] *
      (nx[i] - jags.data$mean.C))
  fec.pred[,i] <- exp(out1$sims.list$alpha[,3] + out1$sims.list$beta[,3] *
      (nx[i] - jags.data$mean.C))
}
quant <- function (x) quantile(x, c(0.025, 0.975))


# Make figure 8.3
library(scales)
library(plotrix)
op <- par(mfrow=c(2,2), mar=c(3,4.5,1.5,2))
plot(NA, ylab=expression('Juvenile survival ('*italic(s)[italic(j)]*')'),
    xlab="", axes=FALSE, ylim=c(0, max(out1$q97.5$phij)),
    xlim=c(min(out1$q2.5$N[1:n]), max(out1$q97.5$N[1:n])))
polygon(x=c(nx, nx.max:nx.min),
    y=c(apply(phij.pred, 2, quant)[1,], apply(phij.pred, 2, quant)[2,length(nx):1]),
    border=NA, col=alpha("salmon2", 0.3))
lines(y=apply(phij.pred, 2, mean), x=nx, lwd=2, col="salmon2")
segments(out1$mean$N[1:n], out1$q2.5$phij, out1$mean$N[1:n], out1$q97.5$phij, col="grey50")
segments(out1$q2.5$N[1:n], out1$mean$phij, out1$q97.5$N[1:n], out1$mean$phij, col="grey50")
points(y=out1$mean$phij, x=out1$mean$N[1:n], pch=16, cex=1.3)
axis(2, las=1)
axis(1)
axis(1, at=c(35, 45, 55, 65, 75), labels=NA, tcl=-0.25)
corner.label("A", font=2, cex=1.5, xoff=2)

plot(NA, ylab=expression('Adult survival ('*italic(s)[italic(a)]*')'),
    xlab="", axes=FALSE, ylim=c(min(out1$q2.5$phia), max(out1$q97.5$phia)),
    xlim=c(min(out1$q2.5$N[1:n]), max(out1$q97.5$N[1:n])))
polygon(x=c(nx, nx.max:nx.min),
    y=c(apply(phia.pred, 2, quant)[1,], apply(phia.pred, 2, quant)[2,length(nx):1]),
    border=NA, col=alpha("salmon2", 0.3))
lines(y=apply(phia.pred, 2, mean), x=nx, lwd=2, col="salmon2")
segments(out1$mean$N[1:n], out1$q2.5$phia, out1$mean$N[1:n], out1$q97.5$phia, col="grey50")
segments(out1$q2.5$N[1:n], out1$mean$phia, out1$q97.5$N[1:n], out1$mean$phia, col="grey50")
points(y=out1$mean$phia, x=out1$mean$N[1:n], pch=16, cex=1.3)
axis(2, las=1)
axis(1)
axis(1, at=c(35, 45, 55, 65, 75), labels=NA, tcl=-0.25)
corner.label("B", font=2, cex=1.5, xoff=2)

par(mar=c(4.5,4.5,0,2))
plot(NA, ylab=expression('Productivity ('*italic(f)*')'),
    xlab=expression('Population size ('*italic(N)*')'), axes=FALSE,
    ylim=c(min(out1$q2.5$f), max(out1$q97.5$f)),
    xlim=c(min(out1$q2.5$N[1:n]), max(out1$q97.5$N[1:n])))
polygon(x=c(nx, nx.max:nx.min),
    y=c(apply(fec.pred, 2, quant)[1,], apply(fec.pred, 2, quant)[2,length(nx):1]),
    border=NA, col=alpha("salmon2", 0.3))
lines(y=apply(fec.pred, 2, mean), x=nx, lwd=2, col="salmon2")
segments(out1$mean$N[1:n], out1$q2.5$f[1:n], out1$mean$N[1:n], out1$q97.5$f[1:n], col="grey50")
segments(out1$q2.5$N[1:n], out1$mean$f[1:n], out1$q97.5$N[1:n], out1$mean$f[1:n], col="grey50")
points(y=out1$mean$f[1:n], x=out1$mean$N[1:n], pch=16, cex=1.3)
axis(2, las=1)
axis(1)
axis(1, at=c(35, 45, 55, 65, 75), labels=NA, tcl=-0.25)
corner.label("C", font=2, cex=1.5, xoff=2)

plot(NA, ylab=expression('Population growth rate ('*lambda*')'),
    xlab=expression('Population size ('*italic(N)*')'), axes=FALSE,
    ylim=c(min(apply(lam, 2, quant)[1,]), max(apply(lam, 2, quant)[2,])),
    xlim=c(min(out1$q2.5$N[1:n]), max(out1$q97.5$N[1:n])))
segments(out1$mean$N[1:n], apply(lam, 2, quant)[1,], out1$mean$N[1:n],
    apply(lam, 2, quant)[2,], col="grey50")
segments(out1$q2.5$N[1:n], apply(lam, 2, mean), out1$q97.5$N[1:n],
    apply(lam, 2, mean), col="grey50")
points(y=apply(lam, 2, mean), x=out1$mean$N[1:n], pch=16, cex=1.3)
axis(2, las=1)
axis(1)
axis(1, at=c(35, 45, 55, 65, 75), labels=NA, tcl=-0.25)
corner.label("D", font=2, cex=1.5, xoff=2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
