# Schaub & Kéry (2022) Integrated Population Models
# Chapter 6 : Benefits of integrated population modeling
# ------------------------------------------------------

# Run time for test script 2 mins, full run 20 mins

# 6.6 Flexibility
# ===============

# 6.6.1 Diversity of data types combined in an IPM (no code)

# 6.6.2 Unequal temporal coverage of data sets – missing values
#       in certain years
# -------------------------------------------------------------

library(IPMbook); library(jagsUI)
data(woodchat66)
marr <- marrayAge(woodchat66$ch, woodchat66$age)

# Bundle data
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], rel.j=rowSums(marr[,,1]),
    rel.a=rowSums(marr[,,2]), C=woodchat66$count, J=woodchat66$J, B=woodchat66$B, n.occasionsT=20,
    n.occasionsC=length(woodchat66$count), n.occasionsCR=dim(marr)[2],
    n.occasionsP=length(woodchat66$J), recap=c(rep(1,3),2,rep(1,5)))

# Write JAGS model file
cat(file="model5.txt", "
model {
  # Priors and linear models
  # Temporal random effects for demographic rates (full time)
  for (t in 1:(n.occasionsT-1)){
    logit(sj[t]) <- mu[1] + eps.sj[t]
    eps.sj[t] ~ dnorm(0, tau.sj)
    logit(sa[t]) <- mu[2] + eps.sa[t]
    eps.sa[t] ~ dnorm(0, tau.sa)
    log(f[t]) <- mu[3] + eps.f[t]
    eps.f[t] ~ dnorm(0, tau.f)
  }

  for (t in 1:(n.occasionsCR-1)){
    p[t] <- mean.p[recap[t]]
  }
  mean.p[1] ~ dunif(0, 1)               # Prior for years with recapture
  mean.p[2] <- 0                        # Fix to zero for year without recaptures

  mu[1] <- logit(mean.sj)
  mu[2] <- logit(mean.sa)
  mu[3] <- log(mean.f)
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)

  sigma.sj ~ dunif(0, 5)
  tau.sj <- pow(sigma.sj, -2)
  sigma.sa ~ dunif(0, 5)
  tau.sa <- pow(sigma.sa, -2)
  sigma.f ~ dunif(0, 5)
  tau.f <- pow(sigma.f, -2)

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time (full time period)
  for (t in 1:(n.occasionsT-1)){
    N[1,t+1] <- f[t] / 2 * sj[t] * (N[1,t] + N[2,t])
    N[2,t+1] <- sa[t] * (N[1,t] + N[2,t])
  }

  # Observation model (for years with counts only)
  for (t in 1:n.occasionsC){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Capture-recapture data (CJS model with multinomial likelihood)
  # For years with CR data only
  for (t in 1:(n.occasionsCR-1)){
    marr.j[t,1:n.occasionsCR] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasionsCR] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasionsCR-1)){
    # Main diagonal
    q[t] <- 1 - p[t]                    # Probability of non-recapture
    pr.j[t,t] <- sj[t+10] * p[t]        # increase index of survival
    pr.a[t,t] <- sa[t+10] * p[t]        # to match with population model
    # Above main diagonal
    for (j in (t+1):(n.occasionsCR-1)){
      pr.j[t,j] <- sj[t+10] * prod(sa[((t+1):j)+10]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(sa[(t:j)+10]) * prod(q[t:(j-1)]) * p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasionsCR-1)){
    pr.j[t,n.occasionsCR] <- 1-sum(pr.j[t,1:(n.occasionsCR-1)])
    pr.a[t,n.occasionsCR] <- 1-sum(pr.a[t,1:(n.occasionsCR-1)])
  }

  # Productivity data (Poisson regression model)
  # For years with productivity data only
  for (t in 1:(n.occasionsP-1)){
    J[t] ~ dpois(f[t+6] * B[t])
  }

  # Derived quantities: total population size (full time period)
  for (t in 1:n.occasionsT){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

# Initial values
inits <- function(){list(mean.p=c(runif(1, 0, 0.5), NA))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "sigma.sj", "sigma.sa", "sigma.f", "sigma",
    "sj", "sa", "f", "N", "Ntot")

# MCMC settings
ni <- 40000; nb <- 10000; nc <- 3; nt <- 30; na <- 2000

# Call JAGS (ART <1 min), check convergence and produce figure
out5 <- jags(jags.data, inits, parameters, "model5.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out5)


# ~~~ code for  Fig. 6.7 ~~~~
av.count <- 1:20
av.count[c(5,6,7,8,15:17,19,20)] <- NA
av.prod <- 1:20
av.prod[c(1:6,16:20)] <- NA
av.cr <- 1:20
av.cr[c(1:10, 15)] <- NA

op <- par(mar=c(4,10,2,1))
plot(y=rep(3, 20), x=av.count, pch=15, type="p", axes=FALSE, ylim=c(0.5, 3.5),
    xlim=c(1,20), xlab="Year", ylab=NA, cex=3)
points(y=rep(2, 20), x=av.prod, pch=15, cex=3)
points(y=rep(1, 20), x=av.cr, pch=15, cex=3)
points(y=1, x=15, pch=0, cex=3)
axis(1, at=1:20, labels=1:20)
axis(2, las=1, at=1:3, labels=c("Capture-recapture", "Productivity",
    "Population counts"))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for Fig. 6.8 ~~~~
op <- par(mar=c(2, 5, 3, 1), las=1, lwd=1, "mfrow")
layout(matrix(1:4, 2, 2, byrow=TRUE), widths=c(1.25, 1.25), heights=c(1, 1.1), TRUE)

u <- col2rgb("grey82")
T <- 20
col.pol <- rgb(u[1], u[2], u[3], alpha=100, maxColorValue=255)

av.count <- rep(16,20)
av.count[c(5,6,7,8,15:17,19,20)] <- 1
plot(out5$mean$Ntot, type="n", ylim=range(c(out5$q2.5$Ntot, out5$q97.5$Ntot)),
    ylab="Population size", xlab=NA, las=1, cex=1.5, axes=FALSE)
axis(1, at=1:T, labels=NA, tcl=-0.25, lwd=1)
axis(1, at=c(5, 10, 15, 20), labels=c(5, 10, 15, 20), tcl=-0.5, lwd=1)
axis(2, las=1, lwd=1)
polygon(c(1:T, T:1), c(out5$q2.5$Ntot, out5$q97.5$Ntot[T:1]), border=NA, col=col.pol)
points(out5$mean$Ntot, type="b", col="black", pch=av.count, lty=1, lwd=1)

x=seq(1.5, 19.5, by=1)
av.prod <- rep(16,20)
av.prod[c(1:6,16:20)] <- 1
plot(y=out5$mean$f, x=x, type="n", ylim=c(1.5, 4.1),
    ylab="Productivity", xlab=NA, axes=FALSE)
axis(1, at=1:T, labels=NA, tcl=-0.25, lwd=1)
axis(1, at=c(5, 10, 15, 20), labels=c(5, 10, 15, 20), tcl=-0.5, lwd=1)
axis(2, las=1, lwd=1)
polygon(c(x, x[(T-1):1]), c(out5$q2.5$f, out5$q97.5$f[(T-1):1]), border=NA, col=col.pol)
points(y=out5$mean$f, x=x, type="b", col="black", pch=av.prod, lty=1, lwd=1)
abline(h=out5$mean$mean.f, col="red", lty=2, lwd=1)

par(mar=c(5, 5, 3, 1), las=1, lwd=1)
av.cr <- rep(16,19)
av.cr[1:10] <- 1
av.cr[15] <- 22
x=seq(1.5, 19.5, by=1)
plot(y=out5$mean$sj, x=x, type="n", ylim=c(0, 0.8),
    ylab="Juvenile survival", xlab="Year", axes=FALSE)
axis(1, at=1:T, labels=NA, tcl=-0.25, lwd=1)
axis(1, at=c(5, 10, 15, 20), labels=c(5, 10, 15, 20), tcl=-0.5, lwd=1)
axis(2, las=1, lwd=1)
polygon(c(x, x[(T-1):1]), c(out5$q2.5$sj, out5$q97.5$sj[(T-1):1]), border=NA, col=col.pol)
points(y=out5$mean$sj, x=x, type="b", col="black", pch=av.cr, lty=1, lwd=1)
abline(h=out5$mean$mean.sj, col="red", lty=2, lwd=1)

x=seq(1.5, 19.5, by=1)
plot(y=out5$mean$sa, x=x, type="n", ylim=c(0, 0.8),
    ylab="Adult survival", xlab="Year", axes=FALSE)
axis(1, at=1:T, labels=NA, tcl=-0.25, lwd=1)
axis(1, at=c(5, 10, 15, 20), labels=c(5, 10, 15, 20), tcl=-0.5, lwd=1)
axis(2, las=1, lwd=1)
polygon(c(x, x[(T-1):1]), c(out5$q2.5$sa, out5$q97.5$sa[(T-1):1]), border=NA, col=col.pol)
points(y=out5$mean$sa, x=x, type="b", col="black", pch=av.cr, lty=1, lwd=1)
abline(h=out5$mean$mean.sa, col="red", lty=2, lwd=1)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 6.6.3 Time points of data collection do not match (no code)
# 6.6.4 Using estimated indices instead of counts for the population-level data (no code)

# 6.6.5 Observation models for the population-level data
# ------------------------------------------------------

# ~~~~ code for the analysis shown in Fig. 6.9 ~~~~
library(IPMbook); library(jagsUI)
data(woodchat64)
marr <- marrayAge(woodchat64$ch, woodchat64$age)

# Bundle data
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat64$count,
    J=woodchat64$J, B=woodchat64$B)

# Write JAGS model code
cat(file="model16.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.occasions-1)){
    logit(sj[t]) <- mu[1] + eps[1,t]
    eps[1,t] ~ dnorm(0, tau[1])
    logit(sa[t]) <- mu[2] + eps[2,t]
    eps[2,t] ~ dnorm(0, tau[2])
    log(f[t]) <- mu[3] + eps[3,t]
    eps[3,t] ~ dnorm(0, tau[3])
    p[t] <- mean.p
  }

  mu[1] <- logit(mean.sj)
  mu[2] <- logit(mean.sa)
  mu[3] <- log(mean.f)

  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)
  mean.p ~ dunif(0, 1)

  # Reparameterize Sigma in terms of SD and rho
  for (i in 1:3){
    tau[i] <- pow(sigma[i], -2)
    sigma[i] ~ dunif(0, 3)
  }

  # Observation error
  sigma.res ~ dunif(0.5, 100)
  tau.res <- pow(sigma.res, -2)

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- f[t] / 2 * sj[t] * (N[1,t] + N[2,t])
    N[2,t+1] <- sa[t] * (N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
    C[t] ~ dnorm(Ntot[t], tau.res)
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
    q[t] <- 1 - p[t]   # Probability of non-recapture
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

  # Productivity data (Poisson regression model)
  for (t in 1:(n.occasions-1)){
    J[t] ~ dpois(f[t] * B[t])
  }
}
")


# Initial values
inits <- function(){list(mean.sa=runif(1, 0.4, 0.6), mean.sj=runif(1, 0.25, 0.35),
    mean.p=runif(1, 0.5, 0.7), mean.f=runif(1, 2.8, 3.3))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "sj", "sa", "f",
    "N", "Ntot", "sigma.res", "sigma")

# MCMC settings
# ni <- 160000; nb <- 110000; nc <- 3; nt <- 10; na <- 2000
ni <- 16000; nb <- 11000; nc <- 3; nt <- 1; na <- 200  # ~~~ for testing

# Call JAGS (ART 6 min)
out16 <- jags(jags.data, inits, parameters, "model16.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# Write JAGS model code
cat(file="model17.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.occasions-1)){
    logit(sj[t]) <- mu[1] + eps[1,t]
    eps[1,t] ~ dnorm(0, tau[1])
    logit(sa[t]) <- mu[2] + eps[2,t]
    eps[2,t] ~ dnorm(0, tau[2])
    log(f[t]) <- mu[3] + eps[3,t]
    eps[3,t] ~ dnorm(0, tau[3])
    p[t] <- mean.p
  }

  mu[1] <- logit(mean.sj)
  mu[2] <- logit(mean.sa)
  mu[3] <- log(mean.f)

  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)
  mean.p ~ dunif(0, 1)

  # Reparameterize Sigma in terms of SD and rho
  for (i in 1:3){
    tau[i] <- pow(sigma[i], -2)
    sigma[i] ~ dunif(0, 3)
  }

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- f[t] / 2 * sj[t] * (N[1,t] + N[2,t])
    N[2,t+1] <- sa[t] * (N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
    C[t] ~ dpois(Ntot[t])
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
    q[t] <- 1 - p[t]   # Probability of non-recapture
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

  # Productivity data (Poisson regression model)
  for (t in 1:(n.occasions-1)){
    J[t] ~ dpois(f[t] * B[t])
  }
}
")

# Initial values
inits <- function(){list(mean.p=runif(1, 0, 0.5))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "sj", "sa", "f",
    "N", "Ntot", "sigma")

# MCMC settings
# ni <- 60000; nb <- 10000; nc <- 3; nt <- 10; na <- 2000
ni <- 6000; nb <- 1000; nc <- 3; nt <- 1; na <- 200  # ~~~ for testing

# Call JAGS from R (ART 5.6 min)
out17 <- jags(jags.data, inits, parameters, "model17.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# Bundle data
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), logC=log(woodchat64$count),
    J= woodchat64$J, B= woodchat64$B)

# Write JAGS model code
cat(file="model18.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.occasions-1)){
    logit(sj[t]) <- mu[1] + eps[1,t]
    eps[1,t] ~ dnorm(0, tau[1])
    logit(sa[t]) <- mu[2] + eps[2,t]
    eps[2,t] ~ dnorm(0, tau[2])
    log(f[t]) <- mu[3] + eps[3,t]
    eps[3,t] ~ dnorm(0, tau[3])
    p[t] <- mean.p
  }

  mu[1] <- logit(mean.sj)
  mu[2] <- logit(mean.sa)
  mu[3] <- log(mean.f)

  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)
  mean.p ~ dunif(0, 1)

  # Reparameterize Sigma in terms of SD and rho
  for (i in 1:3){
    tau[i] <- pow(sigma[i], -2)
    sigma[i] ~ dunif(0, 3)
  }

  # Observation error
  sigma.res ~ dunif(0.001, 10)
  tau.res <- pow(sigma.res, -2)

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- f[t] / 2 * sj[t] * (N[1,t] + N[2,t])
    N[2,t+1] <- sa[t] * (N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
    logC[t] ~ dnorm(log(Ntot[t]), tau.res)
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
    q[t] <- 1 - p[t]   # Probability of non-recapture
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

  # Productivity data (Poisson regression model)
  for (t in 1:(n.occasions-1)){
    J[t] ~ dpois(f[t] * B[t])
  }
}
")

# Initial values
inits <- function(){list(mean.p=runif(1, 0, 0.5))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "sj", "sa", "f",
    "N", "Ntot", "sigma.res", "sigma")

# MCMC settings
# ni <- 60000; nb <- 10000; nc <- 3; nt <- 10; na <- 2000
ni <- 6000; nb <- 1000; nc <- 3; nt <- 1; na <- 200  # ~~~ for testing

# Call JAGS (ART 6 min)
out18 <- jags(jags.data, inits, parameters, "model18.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

save(out16, out17, out18, file = "Data Fig 6.9.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for Fig. 6.9 ~~~~
library(RColorBrewer)
co <- brewer.pal(n=8, name="Blues")[c(8,6,4)]
time <- length(out16$mean$Ntot)
d <- 0.2
limits <- range(c(out16$q2.5$Ntot, out16$q97.5$Ntot, out17$q2.5$Ntot,
    out17$q97.5$Ntot, out18$q2.5$Ntot, out18$q97.5$Ntot))
plot(y=out16$mean$Ntot, x=(1:time)-d, ylim=limits, type="b", pch=16,
    axes=FALSE, ylab="Population size", xlab="Year", col=co[1])
segments((1:time)-d, out16$q2.5$Ntot, (1:time)-d, out16$q97.5$Ntot, col=co[1])
points(y=out17$mean$Ntot, x=1:time, type="b", pch=16, col=co[2])
segments(1:time, out17$q2.5$Ntot, 1:time, out17$q97.5$Ntot, col=co[2])
points(y=out18$mean$Ntot, x=(1:time)+d, type="b", pch=16, col=co[3])
segments((1:time)+d, out18$q2.5$Ntot, (1:time)+d, out18$q97.5$Ntot, col=co[3])
axis(2, las=1)
axis(1, at=1:time, labels=NA, tcl=-0.25)
axis(1, at=seq(1, time, by=2), labels=seq(1, time, by=2), tcl=-0.5)
legend("topleft", pch=rep(16,3), col=co,
    legend=c("Normal", "Poisson", "log-Normal"), bty="n")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for the analysis shown in Fig 6.10 ~~~~
# Bundle data
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat64$count,
    J= woodchat64$J, B= woodchat64$B)

# Write JAGS model code
cat(file="model19.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.occasions-1)){
    logit(sj[t]) <- mu[1] + eps[1,t]
    eps[1,t] ~ dnorm(0, tau[1])
    logit(sa[t]) <- mu[2] + eps[2,t]
    eps[2,t] ~ dnorm(0, tau[2])
    log(f[t]) <- mu[3] + eps[3,t]
    eps[3,t] ~ dnorm(0, tau[3])
    p[t] <- mean.p
  }

  mu[1] <- logit(mean.sj)
  mu[2] <- logit(mean.sa)
  mu[3] <- log(mean.f)

  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)
  mean.p ~ dunif(0, 1)

  # Reparameterize Sigma in terms of SD and rho
  for (i in 1:3){
    tau[i] <- pow(sigma[i], -2)
    sigma[i] ~ dunif(0, 3)
  }

  # Observation error
  sigma.res ~ dunif(0.5, 100)
  tau.res <- pow(sigma.res, -2)

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- f[t] / 2 * sj[t] * N[2,t]
    N[2,t+1] <- sa[t] * (N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
    C[t] ~ dnorm(Ntot[t], tau.res)
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
    q[t] <- 1 - p[t]   # Probability of non-recapture
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

  # Productivity data (Poisson regression model)
  for (t in 1:(n.occasions-1)){
    J[t] ~ dpois(f[t] * B[t])
  }
}
")

# Initial values
inits <- function(){list(mean.p=runif(1, 0, 0.5))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "sj", "sa", "f",
    "N", "Ntot", "sigma.res", "sigma")

# MCMC settings
# ni <- 60000; nb <- 10000; nc <- 3; nt <- 10; na <- 2000
ni <- 6000; nb <- 1000; nc <- 3; nt <- 1; na <- 200  # ~~~ for testing

# Call JAGS from R (ART 5.8 min)
out19 <- jags(jags.data, inits, parameters, "model19.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# Write JAGS model code
cat(file="model20.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.occasions-1)){
    logit(sj[t]) <- mu[1] + eps[1,t]
    eps[1,t] ~ dnorm(0, tau[1])
    logit(sa[t]) <- mu[2] + eps[2,t]
    eps[2,t] ~ dnorm(0, tau[2])
    log(f[t]) <- mu[3] + eps[3,t]
    eps[3,t] ~ dnorm(0, tau[3])
    p[t] <- mean.p
  }

  mu[1] <- logit(mean.sj)
  mu[2] <- logit(mean.sa)
  mu[3] <- log(mean.f)

  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)
  mean.p ~ dunif(0, 1)

  # Reparameterize Sigma in terms of SD and rho
  for (i in 1:3){
    tau[i] <- pow(sigma[i], -2)
    sigma[i] ~ dunif(0, 3)
  }

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- f[t] / 2 * sj[t] * N[2,t]
    N[2,t+1] <- sa[t] * (N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
    C[t] ~ dpois(Ntot[t])
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
    q[t] <- 1 - p[t]   # Probability of non-recapture
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

  # Productivity data (Poisson regression model)
  for (t in 1:(n.occasions-1)){
    J[t] ~ dpois(f[t] * B[t])
  }
}
")

# Initial values
inits <- function(){list(mean.p=runif(1, 0, 0.5))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "sj", "sa", "f",
    "N", "Ntot", "sigma")

# MCMC settings
# ni <- 60000; nb <- 10000; nc <- 3; nt <- 10; na <- 2000
ni <- 6000; nb <- 1000; nc <- 3; nt <- 1; na <- 200  # ~~~ for testing

# Call JAGS from R (ART 5.7 min)
out20 <- jags(jags.data, inits, parameters, "model20.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

# Bundle data
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), logC=log(woodchat64$count),
    J= woodchat64$J, B= woodchat64$B)

# Write JAGS model code
cat(file="model21.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.occasions-1)){
    logit(sj[t]) <- mu[1] + eps[1,t]
    eps[1,t] ~ dnorm(0, tau[1])
    logit(sa[t]) <- mu[2] + eps[2,t]
    eps[2,t] ~ dnorm(0, tau[2])
    log(f[t]) <- mu[3] + eps[3,t]
    eps[3,t] ~ dnorm(0, tau[3])
    p[t] <- mean.p
  }

  mu[1] <- logit(mean.sj)
  mu[2] <- logit(mean.sa)
  mu[3] <- log(mean.f)

  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)
  mean.p ~ dunif(0, 1)

  # Reparameterize Sigma in terms of SD and rho
  for (i in 1:3){
    tau[i] <- pow(sigma[i], -2)
    sigma[i] ~ dunif(0, 3)
  }

  # Observation error
  sigma.res ~ dunif(0.001, 10)
  tau.res <- pow(sigma.res, -2)

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- f[t] / 2 * sj[t] * N[2,t]
    N[2,t+1] <- sa[t] * (N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
    logC[t] ~ dnorm(log(Ntot[t]), tau.res)
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
    q[t] <- 1 - p[t]   # Probability of non-recapture
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

  # Productivity data (Poisson regression model)
  for (t in 1:(n.occasions-1)){
    J[t] ~ dpois(f[t] * B[t])
  }
}
")

# Initial values
inits <- function(){list(mean.p=runif(1, 0, 0.5))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "sj", "sa", "f",
    "N", "Ntot", "sigma.res", "sigma")

# MCMC settings
# ni <- 60000; nb <- 10000; nc <- 3; nt <- 10; na <- 2000
ni <- 6000; nb <- 1000; nc <- 3; nt <- 1; na <- 200  # ~~~ for testing

# Call JAGS from R (ART 5.9 min)
out21 <- jags(jags.data, inits, parameters, "model21.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

count <- woodchat64$count
save(out19, out20, out21, count, file="Data Fig 6.10.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~ code for Fig. 6.10 ~~~~
load("Data Fig 6.9.Rdata")
load("Data Fig 6.10.Rdata")

op <- par(mfrow=c(3,1), las=1, mar=c(3,5,2,1))
boxplot(cbind(out16$sims.list$mean.sj, out17$sims.list$mean.sj,
    out18$sims.list$mean.sj, NA, out19$sims.list$mean.sj,
    out20$sims.list$mean.sj, out21$sims.list$mean.sj),
    outline=FALSE, axes=FALSE, ylab="Juvenile survival", xlab=NA,
    col=c(co, NA, rep("white", 3)), border=c(rep("black", 3), NA, co), boxwex=0.75)
axis(2)
labs <- c("Norm", "Pois", "logNorm", NA, "Norm", "Pois", "logNorm")
axis(1, at=c(1:3, 5:7), labels=NA, tcl=-0.5)
text(x=0.75, y= out20$q97.5$mean.sj, "A", font=2)

boxplot(cbind(out16$sims.list$mean.sa, out17$sims.list$mean.sa,
    out18$sims.list$mean.sa, NA, out19$sims.list$mean.sa,
    out20$sims.list$mean.sa, out21$sims.list$mean.sa),
    outline=FALSE, axes=FALSE, ylab="Adult survival", xlab=NA,
    col=c(co, NA, rep("white", 3)), border=c(rep("black", 3), NA, co), boxwex=0.75)
axis(2)
labs <- c("Norm", "Pois", "logNorm", NA, "Norm", "Pois", "logNorm")
axis(1, at=c(1:3, 5:7), labels=NA, tcl=-0.5)
text(x=0.75, y= out20$q97.5$mean.sa, "B", font=2)

par(mar=c(5,5,2,1))
boxplot(cbind(out16$sims.list$mean.f, out17$sims.list$mean.f,
    out18$sims.list$mean.f, NA, out19$sims.list$mean.f,
    out20$sims.list$mean.f, out21$sims.list$mean.f),
    outline=FALSE, axes=FALSE, ylab="Productivity", xlab=NA,
    col=c(co, NA, rep("white", 3)), border=c(rep("black", 3), NA, co), boxwex=0.75)
axis(2)
labs <- c("Norm", "Pois", "logNorm", NA, "Norm", "Pois", "logNorm")
axis(1, at=c(1:7), labels=labs, tcl=0)
axis(1, at=c(1:3, 5:7), labels=NA, tcl=-0.5)
mtext("Correct IPM", side=1, at=2, line=3)
mtext("Miss-specified IPM", side=1, at=6, line=3)
text(x=0.75, y= out20$q97.5$mean.f, "C", font=2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for Fig 6.11 ~~~~
library(RColorBrewer)
co <- brewer.pal(n=8, name="Blues")[c(8,6,4)]

time <- length(out19$mean$Ntot)
d <- 0.2
limits <- range(c(out19$q2.5$Ntot, out19$q97.5$Ntot, out20$q2.5$Ntot,
    out20$q97.5$Ntot, out21$q2.5$Ntot, out21$q97.5$Ntot))
plot(y= out19$mean$Ntot, x=(1:time)-d, ylim=limits, type="b", pch=16,
    axes=FALSE, ylab="Population size", xlab=NA, col=co[1])
segments((1:time)-d, out19$q2.5$Ntot, (1:time)-d, out19$q97.5$Ntot, col=co[1])
points(y= out20$mean$Ntot, x=1:time, type="b", pch=16, col=co[2])
segments(1:time, out20$q2.5$Ntot, 1:time, out20$q97.5$Ntot, col=co[2])
points(y= out21$mean$Ntot, x=(1:time)+d, type="b", pch=16, col=co[3])
segments((1:time)+d, out21$q2.5$Ntot, (1:time)+d, out21$q97.5$Ntot, col=co[3])
points(woodchat64$count, type = "b", pch = 0)
points(woodchat64$trueN, type = "p", pch = 4, col = "red")
axis(2, las=1)
axis(1, at=1:time, labels=NA, tcl=-0.25)
axis(1, at=seq(1, time, by=2), labels=seq(1, time, by=2), tcl=-0.5)
legend("topright", pch=c(0, 4), col=c("black", "red"),
    legend=c("Observed counts", "True population size"), bty="n", lwd=c(1, NA), title=NA)
legend("top", pch=rep(16,3), col=co,
    legend=c("Normal", "Poisson", "logNormal"), bty="n", title="Estimated population size")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 6.6.6 Informative priors and sequential analyses (no code)
