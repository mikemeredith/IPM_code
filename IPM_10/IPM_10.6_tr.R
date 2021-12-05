# Schaub & Kéry (2022) Integrated Population Models
# Chapter 10 : Population viability analysis
# ------------------------------------------

# Run time for test script 4 mins, full run 36 mins

library(IPMbook) ; library(jagsUI)

# 10.6 PVA of a population with immigration
# =========================================

# Load hoopoe data
library(IPMbook)
data(hoopoe)

# Produce m-array
marr <- marrayAge(hoopoe$ch, hoopoe$age)
marr.j <- marr[,,1]
marr.a <- marr[,,2]

# Bundle data
jags.data <- with(hoopoe$reproAgg, list(marr.j=marr.j, marr.a=marr.a, n.occasions=ncol(marr.j),
  rel.j=rowSums(marr.j), rel.a=rowSums(marr.a), J1=J1, B1=B1, J2=J2, B2=B2, count=hoopoe$count,
  pNinit=dUnif(1, 50), K=10))
str(jags.data)                                    # Remind ourselves of how the data look like

# Write JAGS model file
cat(file="model4.txt", "
model {
  # Priors and linear models
  # For mean and variance parameters
  mean.phij ~ dunif(0, 1)
  mean.phia ~ dunif(0, 1)
  mean.f1 ~ dunif(0, 10)
  mean.f2 ~ dunif(0, 10)
  mean.omega ~ dunif(0, 3)

  mu[1] <- logit(mean.phij)
  mu[2] <- logit(mean.phia)
  mu[3] <- log(mean.omega)
  mu[4] <- log(mean.f1)
  mu[5] <- log(mean.f2)

  mean.p ~ dunif(0, 1)
  mean.pj ~ dunif(0, 1)

  # Linear models for demographic parameters
  for (t in 1:(n.occasions-1+K)){
    logit(phij[t]) <- eps[1,t]
    logit(phia[t]) <- eps[2,t]
    log(omega[t]) <- eps[3,t]
    log(f1[t]) <- eps[4,t]
    log(f2[t]) <- eps[5,t]
    eps[1:5,t] ~ dmnorm.vcov(mu[], sigma2[,])
  }

  for (t in 1:(n.occasions-1)){
    p[t] <- mean.p
    pj[t] <- mean.pj
  }

  # Reparameterize sigma2 in terms of SD and rho
  for (i in 1:5){
    sigma2[i,i] <- sigma[i] * sigma[i]
  }
  for (i in 1:4){
    for (j in (i+1):5){
      sigma2[i,j] <- sigma[i] * sigma[j] * rho[i,j]
      sigma2[j,i] <- sigma2[i,j]
      rho[i,j] ~ dunif(-1, 1)                     # prior for rho
    } #j
  } #i

  # Specify priors for SD
  for (i in 1:5){
    sigma[i] ~ dunif(0, 5)
  }

  sigma.res ~ dunif(0.4, 20)
  tau.res <- pow(sigma.res, -2)

  # Population count data (state-space model)
  # Model for initial stage-spec. population sizes: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)
  N[3,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1+K)){
    N[1,t+1] ~ dpois(f1[t] / 2 * phij[t] * (N[1,t] + N[3,t]) + f2[t] / 2 * phij[t] * N[2,t])
    N[2,t+1] ~ dbin(phia[t], (N[1,t] + N[2,t] + N[3,t]))
    N[3,t+1] ~ dpois((N[1,t] + N[2,t] + N[3,t]) * omega[t])
  }

  # Observation model
  for (t in 1:n.occasions){
    count[t] ~ dnorm(Ntot[t], tau.res)
  }
  for (t in 1:(n.occasions+K)){
    Ntot[t] <- N[1,t] + N[2,t] + N[3,t]           # total population size
  }

  # Productivity data (Poisson regression model)
  for (t in 1:(n.occasions-1)){
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
    qj[t] <- 1-pj[t]
    q[t] <- 1-p[t]
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
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
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
parameters <- c("mean.phij", "mean.phia", "mean.omega", "mean.f1", "mean.f2", "mean.p", "mean.pj",
    "sigma.res", "rho", "sigma2", "phij", "phia", "omega", "f1", "f2", "N", "Ntot", "lambda")

# MCMC settings
# ni <- 110000; nb <- 10000; nc <- 3; nt <- 100; na <- 2000
ni <- 11000; nb <- 1000; nc <- 3; nt <- 10; na <- 200  # ~~~ for testing

# Call JAGS (ART 17 min), check convergence
out4 <- jags(jags.data, inits, parameters, "model4.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out4)


# ~~~~ code for modified model ~~~~
# Write JAGS model file
cat(file="model5.txt", "
model {
  # Priors and linear models
  # For mean and variance parameters
  mean.phij ~ dunif(0, 1)
  mean.phia ~ dunif(0, 1)
  mean.f1 ~ dunif(0, 10)
  mean.f2 ~ dunif(0, 10)
  mean.omega ~ dunif(0, 100)

  mu[1] <- logit(mean.phij)
  mu[2] <- logit(mean.phia)
  mu[3] <- log(mean.omega)
  mu[4] <- log(mean.f1)
  mu[5] <- log(mean.f2)

  mean.p ~ dunif(0, 1)
  mean.pj ~ dunif(0, 1)

  # Linear models for demographic parameters
  for (t in 1:(n.occasions-1+K)){
    logit(phij[t]) <- eps[1,t]
    logit(phia[t]) <- eps[2,t]
    log(omega[t])  <- eps[3,t]
    log(f1[t]) <- eps[4,t]
    log(f2[t]) <- eps[5,t]
    eps[1:5,t] ~ dmnorm.vcov(mu[], sigma2[,])
  }

  for (t in 1:(n.occasions-1)){
    p[t] <- mean.p
    pj[t] <- mean.pj
  }

  # Reparameterize sigma2 in terms of SD and rho
  for (i in 1:5){
    sigma2[i,i] <- sigma[i] * sigma[i]
  }
  for (i in 1:4){
    for (j in (i+1):5){
      sigma2[i,j] <- sigma[i] * sigma[j] * rho[i,j]
      sigma2[j,i] <- sigma2[i,j]
      rho[i,j] ~ dunif(-1,1)
    } #j
  } #i

  # Specify priors for SD and rho
  for (i in 1:5){
    sigma[i] ~ dunif(0, 5)
  }

  sigma.res ~ dunif(0.4, 20)
  tau.res <- pow(sigma.res, -2)

  # Population count data (state-space model)
  # Model for initial stage-spec. population sizes: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)
  N[3,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1+K)){
    N[1,t+1] ~ dpois(f1[t] / 2 * phij[t] * (N[1,t] + N[3,t]) + f2[t] / 2 * phij[t] * N[2,t])
    N[2,t+1] ~ dbin(phia[t], (N[1,t] + N[2,t] + N[3,t]))
    N[3,t+1] ~ dpois(omega[t])
  }

  # Observation model
  for (t in 1:n.occasions){
    count[t] ~ dnorm(Ntot[t], tau.res)
  }
  for (t in 1:(n.occasions+K)){
    Ntot[t] <- N[1,t] + N[2,t] + N[3,t]        # total population size
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
    qj[t] <- 1-pj[t]
    q[t] <- 1-p[t]
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
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
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
parameters <- c("mean.phij", "mean.phia", "mean.omega", "mean.f1", "mean.f2",
    "mean.p", "mean.pj", "sigma.res", "rho", "sigma2", "phij", "phia", "omega",
    "f1", "f2", "N", "Ntot", "lambda")

# MCMC settings
# ni <- 110000; nb <- 10000; nc <- 3; nt <- 100; na <- 2000
ni <- 11000; nb <- 1000; nc <- 3; nt <- 10; na <- 200  # ~~~ for testing

# Call JAGS (ART 18 min), check convergence and summarize posteriors
out5 <- jags(jags.data, inits, parameters, "model5.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out5)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~ save output ~~~
save(out4, out5, file="IPM_10.6_output.RData")
# ~~~~~~~~~~~~~~~~~~~

# ~~~~ Code for Fig 10.8 ~~~~
library(scales)
library(RColorBrewer)
col1 <- brewer.pal(3, 'Blues')
col2 <- brewer.pal(3, 'YlOrRd')
n.years <- length(out5$mean$Ntot)

op <- par(mfrow=c(2, 2), mar=c(4.5, 4.2, 1, 1))
plot(0, 0, ylim=range(c(out4$q2.5$Ntot, out4$q97.5$Ntot, out5$q2.5$Ntot, out5$q97.5$Ntot)),
    xlim=c(0.5, n.years), ylab="Total population size", xlab=NA, las=1,
    col="black", type="l", axes=FALSE)
axis(2, las=1)
axis(1, at=seq(4, n.years, 5), labels=c(2005, 2010, 2015, 2020, 2025))
axis(1, at=1:n.years, labels=NA, tcl=-0.25)
polygon(x=c(1:n.years, n.years:1), y=c(out4$q2.5$Ntot, rev(out4$q97.5$Ntot)),
    col=alpha(col1[2], alpha=0.5), border=NA)
lines(out4$q50$Ntot, col=col1[3], lwd=2)
polygon(x=c(1:n.years, n.years:1), y=c(out5$q2.5$Ntot, rev(out5$q97.5$Ntot)),
    col=alpha(col2[2], alpha=0.5), border=NA)
lines(out5$q50$Ntot, col=col2[3], lwd=2)
legend("topleft", col=c(col1[3], col2[3]), lwd=rep(2, 2),
    legend=c("rate", "number"), bty="n", title="Parameterization of immigration")

plot(0, 0, ylim=range(c(out4$q2.5$N[1,], out4$q97.5$N[1,], out5$q2.5$N[1,], out5$q97.5$N[1,])),
    xlim=c(0.5, n.years), ylab="Local recruits", xlab=NA, las=1,
    col="black", type="l", axes=FALSE)
axis(2, las=1)
axis(1, at=seq(4, n.years, 5), labels=c(2005, 2010, 2015, 2020, 2025))
axis(1, at=1:n.years, labels=NA, tcl=-0.25)
polygon(x=c(1:n.years, n.years:1), y=c(out4$q2.5$N[1,], rev(out4$q97.5$N[1,])),
    col=alpha(col1[2], alpha=0.5), border=NA)
lines(out4$q50$N[1,], col=col1[3], lwd=2)
polygon(x=c(1:n.years, n.years:1), y=c(out5$q2.5$N[1,], rev(out5$q97.5$N[1,])),
    col=alpha(col2[2], alpha=0.5), border=NA)
lines(out5$q50$N[1,], col=col2[3], lwd=2)

plot(0, 0, ylim=range(c(out4$q2.5$N[2,], out4$q97.5$N[2,], out5$q2.5$N[2,], out5$q97.5$N[2,])),
    xlim=c(0.5, n.years), ylab="Surviving adults", xlab=NA, las=1,
    col="black", type="l", axes=FALSE)
axis(2, las=1)
axis(1, at=seq(4, n.years, 5), labels=c(2005, 2010, 2015, 2020, 2025))
axis(1, at=1:n.years, labels=NA, tcl=-0.25)
polygon(x=c(1:n.years, n.years:1), y=c(out4$q2.5$N[2,], rev(out4$q97.5$N[2,])),
    col=alpha(col1[2], alpha=0.5), border=NA)
lines(out4$q50$N[2,], col=col1[3], lwd=2)
polygon(x=c(1:n.years, n.years:1), y=c(out5$q2.5$N[2,], rev(out5$q97.5$N[2,])),
    col=alpha(col2[2], alpha=0.5), border=NA)
lines(out5$q50$N[2,], col=col2[3], lwd=2)

plot(0, 0, ylim=range(c(out4$q2.5$N[3,], out4$q97.5$N[3,], out5$q2.5$N[3,], out5$q97.5$N[3,])),
    xlim=c(0.5, n.years), ylab="Number of immigrants", xlab=NA, las=1,
    col="black", type="l", axes=FALSE)
axis(2, las=1)
axis(1, at=seq(4, n.years, 5), labels=c(2005, 2010, 2015, 2020, 2025))
axis(1, at=1:n.years, labels=NA, tcl=-0.25)
polygon(x=c(1:n.years, n.years:1), y=c(out4$q2.5$N[3,], rev(out4$q97.5$N[3,])),
    col=alpha(col1[2], alpha=0.5), border=NA)
lines(out4$q50$N[3,], col=col1[3], lwd=2)
polygon(x=c(1:n.years, n.years:1), y=c(out5$q2.5$N[3,], rev(out5$q97.5$N[3,])),
    col=alpha(col2[2], alpha=0.5), border=NA)
lines(out5$q50$N[3,], col=col2[3], lwd=2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
