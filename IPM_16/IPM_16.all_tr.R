# Schaub & Kéry (2022) Integrated Population Models
# Chapter 16 : Barn swallow
# -------------------------

# Run time for test script 5 mins, full run 40 mins

# 16.4 Component data likelihoods
# =============================================

library(IPMbook); library(jagsUI)
data(swallow)
str(swallow)
# List of 5
# $ marr.j      : num [1:9, 1:6, 1:7] 0 0 0 0 1 NA NA NA NA 0 ...
# ..- attr(*, "dimnames")=List of 3
# .. ..$ : chr [1:9] "Buus" "Baulmes" "Tavannes" "Riviera" ...
# .. ..$ : chr [1:6] "1997" "1998" "1999" "2000" ...
# .. ..$ : chr [1:7] "1998" "1999" "2000" "2001" ...
# $ marr.a      : num [1:9, 1:6, 1:7] 19 7 9 8 9 NA NA NA NA 0 ...
# ..- attr(*, "dimnames")=List of 3
# .. ..$ : chr [1:9] "Buus" "Baulmes" "Tavannes" "Riviera" ...
# .. ..$ : chr [1:6] "1997" "1998" "1999" "2000" ...
# .. ..$ : chr [1:7] "1998" "1999" "2000" "2001" ...
# $ counts      : num [1:9, 1:7] 75 NA 43 NA 39 NA NA NA NA 77 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:9] "Buus" "Baulmes" "Tavannes" "Riviera" ...
# .. ..$ : chr [1:7] "1997" "1998" "1999" "2000" ...
# $ productivity: num [1:2, 1:9, 1:7, 1:3] 359 208 233 79 183 93 92 48 ...
# ..- attr(*, "dimnames")=List of 4
# .. ..$ : chr [1:2] "First" "Second"
# .. ..$ : chr [1:9] "Buus" "Baulmes" "Tavannes" "Riviera" ...
# .. ..$ : chr [1:7] "1997" "1998" "1999" "2000" ...
# .. ..$ : chr [1:3] "Eggs" "Fledglings" "Broods"
# $ second      : num [1:9, 1:7, 1:2] 77 52 39 20 17 0 0 0 0 80 ...
# ..- attr(*, "dimnames")=List of 3
# .. ..$ : chr [1:9] "Buus" "Baulmes" "Tavannes" "Riviera" ...
# .. ..$ : chr [1:7] "1997" "1998" "1999" "2000" ...
# .. ..$ : chr [1:2] "First broods" "Second broods"


# 16.4.1 Population count data (no code)
# 16.4.2 Productivity data (no code)
# 16.4.3 Capture-recapture data (no code)

# 16.5 The integrated population model
# ====================================

# Define the starting year of each site
start <- c(1, 1, 1, 1, 1, 3, 3, 3, 4)

# Compute the total number of released individuals
rel.j <- apply(swallow$marr.j, c(1,2), sum, na.rm=TRUE)
rel.a <- apply(swallow$marr.a, c(1,2), sum, na.rm=TRUE)

# Define ranges of discrete uniform priors for initial stage-structured population sizes per site
pinit <- dUnif(lower=rep(1,9), upper=apply(swallow$counts, 1, max, na.rm=TRUE))

# Bundle data and produce data overview
jags.data <- with(swallow, list(marr.j=marr.j, rel.j=rel.j, marr.a=marr.a, rel.a=rel.a, y=counts,
    c=productivity[,,,'Eggs'], b=productivity[,,,'Broods'], f=productivity[,,,'Fledglings'],
    v=second[,,'First broods'], w=second[,,'Second broods'], pinit=pinit, nsites=9,
    nyears=7, start=start))
str(jags.data)
# List of 14
# $ marr.j: num [1:9, 1:6, 1:7] 0 0 0 0 1 NA NA NA NA 0 ...
# $ rel.j : num [1:9, 1:6] 229 161 128 66 92 0 0 0 0 269 ...
# $ marr.a: num [1:9, 1:6, 1:7] 19 7 9 8 9 NA NA NA NA 0 ...
# $ rel.a : num [1:9, 1:6] 57 44 31 19 29 0 0 0 0 62 ...
# $ y     : num [1:9, 1:7] 75 NA 43 NA 39 NA NA NA NA 77 ...
# $ c     : num [1:2, 1:9, 1:7] 359 208 233 79 183 93 92 48 71 4 ...
# $ b     : num [1:2, 1:9, 1:7] 77 47 52 18 39 21 20 12 17 1 ...
# $ f     : num [1:2, 1:9, 1:7] 281 179 178 61 148 78 66 37 60 4 ...
# $ v     : num [1:9, 1:7] 77 52 39 20 17 0 0 0 0 80 ...
# $ w     : num [1:9, 1:7] 47 18 21 12 1 0 0 0 0 55 ...
# $ pinit : num [1:9, 1:115] 0.0116 0.013 0.0227 0.0161 0.0256 ...
# $ nsites: num 9
# $ nyears: num 7
# $ start : num [1:9] 1 1 1 1 1 3 3 3 4

# Write JAGS model file
cat(file = "model1.txt", "
model {
  # Priors and linear models
  # Linear models for the eight demographic parameters
  for (s in 1:nsites){
    for (t in start[s]:(nyears-1)){
      logit(phij[s,t]) <- mu.phij[s] + eps.phij[t] + eta.phij[s,t]
      logit(phia[s,t]) <- mu.phia[s] + eps.phia[t] + eta.phia[s,t]
      log(omega[s,t]) <- mu.omega[s] + eps.omega[t] + eta.omega[s,t]
      pj[s,t] ~ dunif(0, 1)
      pa[s,t] ~ dunif(0, 1)
    } # t
    for (t in start[s]:nyears){
      log(rho[1,s,t]) <- mu.rho[1,s] + eps.rho[1,t] + eta.rho[1,s,t]
      log(rho[2,s,t]) <- mu.rho[2,s] + eps.rho[2,t] + eta.rho[2,s,t]
      logit(zeta[1,s,t]) <- mu.zeta[1,s] + eps.zeta[1,t] + eta.zeta[1,s,t]
      logit(zeta[2,s,t]) <- mu.zeta[2,s] + eps.zeta[2,t] + eta.zeta[2,s,t]
      logit(kappa[s,t]) <- mu.kappa[s] + eps.kappa[t] + eta.kappa[s,t]
    } #t
  } #s

  # Priors for the means of the eight demographic parameters
  for (s in 1:nsites){
    mu.phij[s] ~ dnorm(0, 0.001)
    overall.mean[1,s] <- ilogit(mu.phij[s])
    mu.phia[s] ~ dnorm(0, 0.001)
    overall.mean[2,s] <- ilogit(mu.phia[s])
    mu.omega[s] ~ dnorm(0, 0.001)
    overall.mean[3,s] <- exp(mu.omega[s])
    mu.rho[1,s] ~ dnorm(0, 0.001)
    overall.mean[4,s] <- exp(mu.rho[1,s])
    mu.rho[2,s] ~ dnorm(0, 0.001)
    overall.mean[5,s] <- exp(mu.rho[2,s])
    mu.zeta[1,s] ~ dnorm(0, 0.001)
    overall.mean[6,s] <- ilogit(mu.zeta[1,s])
    mu.zeta[2,s] ~ dnorm(0, 0.001)
    overall.mean[7,s] <- ilogit(mu.zeta[2,s])
    mu.kappa[s] ~ dnorm(0, 0.001)
    overall.mean[8,s] <- ilogit(mu.kappa[s])
  }

  # Priors for temporal random effects (common across sites)
  for (t in 1:(nyears-1)){
    eps.phij[t] ~ dnorm(0, tau.t[1])
    eps.phia[t] ~ dnorm(0, tau.t[2])
    eps.omega[t] ~ dnorm(0, tau.t[3])
  }
  for (t in 1:nyears){
    eps.rho[1,t] ~ dnorm(0, tau.t[4])
    eps.rho[2,t] ~ dnorm(0, tau.t[5])
    eps.zeta[1,t] ~ dnorm(0, tau.t[6])
    eps.zeta[2,t] ~ dnorm(0, tau.t[7])
    eps.kappa[t] ~ dnorm(0, tau.t[8])
  }

  # Priors for within-site temporal random effects
  for (s in 1:nsites){
    for (t in start[s]:(nyears-1)){
      eta.phij[s,t] ~ dnorm(0, tau.st[1])
      eta.phia[s,t] ~ dnorm(0, tau.st[2])
      eta.omega[s,t] ~ dnorm(0, tau.st[3])
    } #t
    for (t in start[s]:nyears){
      eta.rho[1,s,t] ~ dnorm(0, tau.st[4])
      eta.rho[2,s,t] ~ dnorm(0, tau.st[5])
      eta.zeta[1,s,t] ~ dnorm(0, tau.st[6])
      eta.zeta[2,s,t] ~ dnorm(0, tau.st[7])
      eta.kappa[s,t] ~ dnorm(0, tau.st[8])
    } #t
  } #s

  # Priors for precisions of the eight demographic parameters
  for (i in 1:8){
    tau.t[i] ~ dgamma(0.001, 0.001)
    tau2[i] <- 1 / tau.t[i]
    tau.st[i] ~ dgamma(0.001, 0.001)
    ypsilon2[i] <- 1 / tau.st[i]
    # Intraclass correlation coefficients (ICC)
    ICC[i] <- tau2[i] / (tau2[i] + ypsilon2[i])
  }
  # Observation error for state-space model
  for (s in 1:nsites){
    tau.y[s] ~ dgamma(0.001, 0.001)
    sigma2[s] <- 1 / tau.y[s]
  }

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  for (s in 1:nsites){
    R[s,start[s]] ~ dcat(pinit[s,]) # Local recruits
    S[s,start[s]] ~ dcat(pinit[s,]) # Surviving adults
    I[s,start[s]] ~ dcat(pinit[s,]) # Immigrants
  }

  # Process model over time: our model of population dynamics
  for (s in 1:nsites){
    for (t in (start[s]+1):nyears){
      R[s,t] ~ dpois((rho[1,s,t-1] * zeta[1,s,t-1] + kappa[s,t-1] * rho[2,s,t-1] * zeta[2,s,t-1]) *
          0.5 * phij[s,t-1] * N[s,t-1])
      S[s,t] ~ dbin(phia[s,t-1], N[s,t-1])
      I[s,t] ~ dpois(N[s,t-1] * omega[s,t-1])
    } #t
  } #s

  # Observation model
  for (s in 1:nsites){
    for(t in start[s]:nyears){
      N[s,t] <- R[s,t] + S[s,t] + I[s,t]
      logN[s,t] <- log(N[s,t])
      y[s,t] ~ dlnorm(logN[s,t], tau.y[s])
    } #t
  } #s

  # Productivity data (Poisson and binomial models)
  for (s in 1:nsites){
    for (t in start[s]:nyears){
      # Clutch size and fledging success of first brood
      c[1,s,t] ~ dpois(rho[1,s,t] * b[1,s,t])
      f[1,s,t] ~ dbin(zeta[1,s,t], c[1,s,t])
      # Clutch size and fledging success of second brood
      c[2,s,t] ~ dpois(rho[2,s,t] * b[2,s,t])
      f[2,s,t] ~ dbin(zeta[2,s,t], c[2,s,t])
      # Probability to conduct a second brood
      w[s,t] ~ dbin(kappa[s,t], v[s,t])
    } #t
  } #s

  # Capture-recapture data (CJS model with multinomial likelihood)
  for (s in 1:nsites){
    # Define the multinomial likelihood
    for (t in start[s]:(nyears-1)){
      marr.j[s,t,start[s]:nyears] ~ dmulti(pr.j[s,t,start[s]:nyears], rel.j[s,t])
      marr.a[s,t,start[s]:nyears] ~ dmulti(pr.a[s,t,start[s]:nyears], rel.a[s,t])
    } #t
    # Define the cell probabilities of the m-arrays
    for (t in start[s]:(nyears-1)){
      # Main diagonal
      qj[s,t] <- 1-pj[s,t]
      qa[s,t] <- 1-pa[s,t]
      pr.j[s,t,t] <- phij[s,t] * pj[s,t]
      pr.a[s,t,t] <- phia[s,t] * pa[s,t]
      # Above main diagonal
      for (j in (t+1):(nyears-1)){
        pr.j[s,t,j] <- phij[s,t] * prod(phia[s,(t+1):j]) * qj[s,t] * prod(qa[s,t:(j-1)]) *
            pa[s,j] / qa[s,t]
        pr.a[s,t,j] <- prod(phia[s,t:j]) * prod(qa[s,t:(j-1)]) * pa[s,j]
      } #j
      # Below main diagonal
      for (j in 1:(t-1)){
        pr.j[s,t,j] <- 0
        pr.a[s,t,j] <- 0
      } #j
      # Last column: probability of non-recapture
      pr.j[s,t,nyears] <- 1-sum(pr.j[s,t,start[s]:(nyears-1)])
      pr.a[s,t,nyears] <- 1-sum(pr.a[s,t,start[s]:(nyears-1)])
    } #t
  } #s
}
")

# Initial values
inits <- function(){list(mu.phia=runif(9, -1, 0))}

# Parameters monitored
parameters <- c("phij", "phia", "omega", "rho", "zeta", "kappa", "pj", "pa", "R", "S", "I", "N",
    "overall.mean", "eps.phij", "eps.phia", "eps.omega", "eps.rho", "eps.zeta", "eps.kappa", "eta.phij",
    "eta.phia", "eta.omega", "eta.rho", "eta.zeta", "eta.kappa", "tau2", "ypsilon2", "ICC", "sigma2")

# MCMC settings
# ni <- 50000; nb <- 10000; nc <- 3; nt <- 40; na <- 5000
ni <- 5000; nt <- 4; nb <- 1000; nc <- 3; na <- 500  # ~~~ for testing

# Call JAGS from R (ART 12 min) and check convergence
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1)


# ICC for the population growth rate
# ''''''''''''''''''''''''''''''''''

# Bundle data
# Calculate population growth rate
r <- log(out1$sims.list$N[,,2:7]) - log(out1$sims.list$N[,,1:6])
jags.data <- list(r=r, start=start, nsites=9, nyears=6, nsamp=dim(r)[1])

cat(file = "model2.txt", "
model {
  # Priors and linear models
  # Fixed site effects
  for (s in 1:nsites){
    mu[s] ~ dnorm(0, 0.001)
  }

  # Temporal random effects
  for (t in 1:nyears){
    epsilon[t] ~ dnorm(0, tau.t)
  }

  # Priors for precisions
  tau.t ~ dgamma(0.001, 0.001)
  tau.ts ~ dgamma(0.001, 0.001)                   # Residual, site-time variation

  # Likelihood
  for (i in 1:nsamp){
    for (s in 1:nsites){
      for (t in start[s]:nyears){
        r[i,s,t] ~ dnorm(mu[s] + epsilon[t], tau.ts)
      } #t
    } #s
  } #i

  # Derived quantities
  # Variances
  tau2 <- 1 / tau.t
  ypsilon2 <- 1 / tau.ts

  # Intraclass correlation coefficient (ICC)
  ICC <- tau2 / (tau2 + ypsilon2)
}
")

# Initial values
inits <- function(){list(tau.t=0.01)}

# Parameters monitored
parameters <- c("ICC", "tau2", "ypsilon2", "mu")

# MCMC settings
# ni <- 10000; nb <- 5000; nc <- 3; nt <- 1; na <- 2000
ni <- 1000; nt <- 1; nb <- 500; nc <- 3; na <- 200  # ~~~ for testing

# Call JAGS from R (ART 54 min) and check convergence
out2 <- jags(jags.data, inits, parameters, "model2.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out2)


# 16.6 Results
# ============

print(out1, 3)

                  # mean     sd     2.5%      50%    97.5% overlap0     f  Rhat n.eff
# phij[1,1]        0.004  0.003    0.000    0.003    0.011    FALSE 1.000 1.003   759
# phij[2,1]        0.005  0.004    0.000    0.004    0.014    FALSE 1.000 1.003  3000
# phij[3,1]        0.000  0.000    0.000    0.000    0.001    FALSE 1.000 1.089  1624
# phij[4,1]        0.006  0.005    0.000    0.005    0.018    FALSE 1.000 1.001  2262
# phij[5,1]        0.011  0.008    0.002    0.009    0.031    FALSE 1.000 1.002  3000
# phij[1,2]        0.005  0.003    0.001    0.004    0.013    FALSE 1.000 1.004   913
# phij[2,2]        0.007  0.006    0.001    0.006    0.023    FALSE 1.000 1.001  3000
# phij[3,2]        0.000  0.000    0.000    0.000    0.001    FALSE 1.000 1.029  2816
# phij[4,2]        0.009  0.007    0.001    0.007    0.026    FALSE 1.000 1.002  1658
# phij[5,2]        0.012  0.007    0.003    0.010    0.031    FALSE 1.000 1.000  3000
# phij[1,3]        0.005  0.003    0.001    0.004    0.012    FALSE 1.000 1.002  1109
# phij[2,3]        0.005  0.004    0.000    0.004    0.013    FALSE 1.000 1.001  1383
# phij[3,3]        0.000  0.000    0.000    0.000    0.001    FALSE 1.000 1.026  3000
# phij[4,3]        0.006  0.005    0.001    0.005    0.017    FALSE 1.000 1.002  3000
# [... output truncated ...]

print(out2, 3)

             # mean      sd      2.5%       50%     97.5% overlap0     f  Rhat n.eff
# ICC         0.320   0.147     0.120     0.288     0.685    FALSE 1.000 1.000 15000
# tau2        0.011   0.017     0.003     0.008     0.041    FALSE 1.000 1.056 13423
# ypsilon2    0.019   0.000     0.019     0.019     0.019    FALSE 1.000 1.000     1
# mu[1]      -0.099   0.044    -0.184    -0.099    -0.014    FALSE 0.984 1.001 15000
# mu[2]      -0.138   0.044    -0.224    -0.138    -0.053    FALSE 0.996 1.001 15000
# mu[3]      -0.045   0.044    -0.130    -0.045     0.041     TRUE 0.880 1.001 15000
# mu[4]      -0.173   0.044    -0.258    -0.173    -0.088    FALSE 0.998 1.001 15000
# mu[5]      -0.128   0.044    -0.214    -0.128    -0.043    FALSE 0.995 1.001 15000
# mu[6]      -0.247   0.044    -0.332    -0.247    -0.162    FALSE 1.000 1.001 15000
# mu[7]      -0.120   0.044    -0.205    -0.120    -0.034    FALSE 0.992 1.001 15000
# mu[8]      -0.078   0.044    -0.163    -0.078     0.007     TRUE 0.968 1.001 15000
# mu[9]      -0.217   0.044    -0.302    -0.217    -0.132    FALSE 0.999 1.001 15000


# ~~~~ code for Fig. 16.4 ~~~~
library(scales)
cl <- viridis_pal(option='E')(9)
qu <- function(x) quantile(x, c(0.025, 0.975), na.rm=TRUE)
ptc <- 1.2
op <- par(las=1, mfrow=c(2, 1), mar=c(2,4,2,1))

# Population size
d <- seq(-0.3, 0.3, length.out=9)
xmean <- seq(1, 7, by=1)
plot(x=xmean+d[1], y=out1$mean$N[1,], type='b', pch=16, ylab='Population size',
    ylim=c(0, 250), xlab=NA, xlim=c(0.7,7.3), axes=FALSE, col=cl[1], cex=ptc)
segments(xmean+d[1], out1$q2.5$N[1,], xmean+d[1], out1$q97.5$N[1,], col=cl[1])
for (s in 2:9){
  points(x=xmean+d[s], y=out1$mean$N[s,], type='b', pch=16, col=cl[s], cex=ptc)
  segments(xmean+d[s], out1$q2.5$N[s,], xmean+d[s], out1$q97.5$N[s,], col=cl[s])
}
axis(2)
axis(1, at=1:7, labels=NA)
legend(x=4.2, y=250, pch=rep(16,4), col=cl[1:4],
    legend=c('Buus', 'Baulmes', 'Tavannes', 'Riviera'), bty='n', pt.cex=ptc)
legend(x=6, y=250, pch=rep(16,5), col=cl[5:9],
    legend=c('Sargans', 'Wauwil', 'Vaulruz', 'Enhaut', 'Dompierre'), bty='n', pt.cex=ptc)

# Calculate population growth rate
lam <- out1$sims.list$N[,,2:7] / out1$sims.list$N[,,1:6]

xmean <- seq(1.5, 6.5, by=1)
par(mar=c(4,4,0,1))
plot(x=xmean+d[1], y=apply(lam, c(2,3), median)[1,], type='b', pch=16,
    ylab='Population growth rate', ylim=c(0.4, 1.4), xlab=NA, xlim=c(0.7,7.3),
    axes=FALSE, col=cl[1], cex=ptc)
segments(xmean+d[1], apply(lam, c(2,3), qu)[1,1,], xmean+d[1],
    apply(lam, c(2,3), qu)[2,1,], col=cl[1])
for (s in 2:9){
  points(x=xmean+d[s], y=apply(lam, c(2,3), mean)[s,], type='b', pch=16, col=cl[s], cex=ptc)
  segments(xmean+d[s], apply(lam, c(2,3), qu)[1,s,], xmean+d[s],
      apply(lam, c(2,3), qu)[2,s,], col=cl[s])
}
axis(2)
axis(1, at=1:7, labels=1997:2003)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 16.5 ~~~~
library(scales)
cl <- viridis_pal(option='E')(9)

op <- par(las=1, mfrow=c(4, 2), mar=c(1, 4.2, 2, 1))
d <- seq(-0.3, 0.3, length.out=9)
ptc <- 1.2

# Juvenile survival
xmean <- seq(1.5, 6.5, by=1)
plot(x=xmean+d[1], y=out1$mean$phij[1,], type='b', pch=16,
    ylab=expression(paste('Juvenile survival (', phi[italic(j)], ')')),
    ylim=c(0, 0.08), xlab=NA, xlim=c(0.7,7.3), axes=FALSE, col=cl[1], cex=ptc)
segments(xmean+d[1], out1$q2.5$phij[1,], xmean+d[1], out1$q97.5$phij[1,], col=cl[1])
for (s in 2:9){
  points(x=xmean+d[s], y=out1$mean$phij[s,], type='b', pch=16, col=cl[s], cex=ptc)
  segments(xmean+d[s], out1$q2.5$phij[s,], xmean+d[s], out1$q97.5$phij[s,], col=cl[s])
}
axis(2)
axis(1, at=1:7, labels=NA)
legend(x=0.7, y=0.085, pch=rep(16,4), col=cl[1:4],
    legend=c('Buus', 'Baulmes', 'Tavannes', 'Riviera'), bty='n', pt.cex=ptc)
legend(x=2.7, y=0.085, pch=rep(16,5), col=cl[5:9],
    legend=c('Sargans', 'Wauwil', 'Vaulruz', 'Enhaut', 'Dompierre'), bty='n', pt.cex=ptc)

# Adult survival
xmean <- seq(1.5, 6.5, by=1)
plot(x=xmean+d[1], y=out1$mean$phia[1,], type='b', pch=16,
    ylab=expression(paste('Adult survival (', phi[italic(a)], ')')),
    ylim=c(0.15, 0.75), xlab=NA, xlim=c(0.7,7.3), axes=FALSE, col=cl[1], cex=ptc)
segments(xmean+d[1], out1$q2.5$phia[1,], xmean+d[1], out1$q97.5$phia[1,], col=cl[1])
for (s in 2:9){
  points(x=xmean+d[s], y=out1$mean$phia[s,], type='b', pch=16, col=cl[s], cex=ptc)
  segments(xmean+d[s], out1$q2.5$phia[s,], xmean+d[s], out1$q97.5$phia[s,], col=cl[s])
}
axis(2)
axis(1, at=1:7, labels=NA)

# Clutch size 1.brood
xmean <- seq(1, 7, by=1)
plot(x=xmean+d[1], y=out1$mean$rho[1,1,], type='b', pch=16,
    ylab=expression(paste('Clutch size 1. brood (', rho[italic(1)], ')')),
    ylim=c(3.5, 5.5), xlab=NA, xlim=c(0.7,7.3), axes=FALSE, col=cl[1], cex=ptc)
segments(xmean+d[1], out1$q2.5$rho[1,1,], xmean+d[1], out1$q97.5$rho[1,1,], col=cl[1])
for (s in 2:9){
  points(x=xmean+d[s], y=out1$mean$rho[1,s,], type='b', pch=16, col=cl[s], cex=ptc)
  segments(xmean+d[s], out1$q2.5$rho[1,s,], xmean+d[s], out1$q97.5$rho[1,s,], col=cl[s])
}
axis(2)
axis(1, at=1:7, labels=NA)

# Clutch size 2.brood
xmean <- seq(1, 7, by=1)
plot(x=xmean+d[1], y=out1$mean$rho[2,1,], type='b', pch=16,
    ylab=expression(paste('Clutch size 2. brood (', rho[italic(2)], ')')),
    ylim=c(3.5, 5.5), xlab=NA, xlim=c(0.7,7.3), axes=FALSE, col=cl[1], cex=ptc)
segments(xmean+d[1], out1$q2.5$rho[2,1,], xmean+d[1], out1$q97.5$rho[2,1,], col=cl[1])
for (s in 2:9){
  points(x=xmean+d[s], y=out1$mean$rho[2,s,], type='b', pch=16, col=cl[s], cex=ptc)
  segments(xmean+d[s], out1$q2.5$rho[2,s,], xmean+d[s], out1$q97.5$rho[2,s,], col=cl[s])
}
axis(2)
axis(1, at=1:7, labels=NA)

# Fledging success 1.brood
xmean <- seq(1, 7, by=1)
plot(x=xmean+d[1], y=out1$mean$zeta[1,1,], type='b', pch=16,
    ylab=expression(paste('Fledging success 1. brood (', zeta[italic(1)], ')')),
    ylim=c(0.4, 1), xlab=NA, xlim=c(0.7,7.3), axes=FALSE, col=cl[1], cex=ptc)
segments(xmean+d[1], out1$q2.5$zeta[1,1,], xmean+d[1], out1$q97.5$zeta[1,1,], col=cl[1])
for (s in 2:9){
  points(x=xmean+d[s], y=out1$mean$zeta[1,s,], type='b', pch=16, col=cl[s], cex=ptc)
  segments(xmean+d[s], out1$q2.5$zeta[1,s,], xmean+d[s], out1$q97.5$zeta[1,s,], col=cl[s])
}
axis(2)
axis(1, at=1:7, labels=NA)

# Fledging success 2.brood
xmean <- seq(1, 7, by=1)
plot(x=xmean+d[1], y=out1$mean$zeta[2,1,], type='b', pch=16,
    ylab=expression(paste('Fledging success 2. brood (', zeta[italic(2)], ')')),
    ylim=c(0.4, 1), xlab=NA, xlim=c(0.7,7.3), axes=FALSE, col=cl[1], cex=ptc)
segments(xmean+d[1], out1$q2.5$zeta[2,1,], xmean+d[1], out1$q97.5$zeta[2,1,], col=cl[1])
for (s in 2:9){
  points(x=xmean+d[s], y=out1$mean$zeta[2,s,], type='b', pch=16, col=cl[s], cex=ptc)
  segments(xmean+d[s], out1$q2.5$zeta[2,s,], xmean+d[s], out1$q97.5$zeta[2,s,], col=cl[s])
}
axis(2)
axis(1, at=1:7, labels=NA)

par(mar=c(3, 4.2, 2, 1))
# Immigration rate
xmean <- seq(1.5, 6.5, by=1)
plot(x=xmean+d[1], y=out1$mean$omega[1,], type='b', pch=16,
    ylab=expression(paste('Immigration rate (', omega, ')')),
    ylim=c(0.15, 0.9), xlab=NA, xlim=c(0.7,7.3), axes=FALSE, col=cl[1], cex=ptc)
segments(xmean+d[1], out1$q2.5$omega[1,], xmean+d[1], out1$q97.5$omega[1,], col=cl[1])
for (s in 2:9){
  points(x=xmean+d[s], y=out1$mean$omega[s,], type='b', pch=16, col=cl[s], cex=ptc)
  segments(xmean+d[s], out1$q2.5$omega[s,], xmean+d[s], out1$q97.5$omega[s,], col=cl[s])
}
axis(2)
axis(1, at=1:7, labels=1997:2003)

# Probability of double brooding
xmean <- seq(1, 7, by=1)
plot(x=xmean+d[1], y=out1$mean$kappa[1,], type='b', pch=16,
    ylab=expression(paste('Double brooding (', kappa, ')')),
    ylim=c(0.3, 0.85), xlab=NA, xlim=c(0.7,7.3), axes=FALSE, col=cl[1], cex=ptc)
segments(xmean+d[1], out1$q2.5$kappa[1,], xmean+d[1], out1$q97.5$kappa[1,], col=cl[1])
for (s in 2:9){
  points(x=xmean+d[s], y=out1$mean$kappa[s,], type='b', pch=16, col=cl[s], cex=ptc)
  segments(xmean+d[s], out1$q2.5$kappa[s,], xmean+d[s], out1$q97.5$kappa[s,], col=cl[s])
}
axis(2)
axis(1, at=1:7, labels=1997:2003)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 16.6 ~~~~
library(denstrip)

order <- c(9, 8, 3, 7, 6, 5, 4, 2)
op <- par(mar=c(4.5, 11, 1, 1))
plot(0, ylim=c(0.8, 9.2), xlim = c(0, 1), axes=FALSE, pch=NA,
    xlab="Intra-class correlation", ylab=NA)
median <- out1$q50$ICC
median[1] <- NA
for (i in 1:8){
  denstrip(out1$sims.list$ICC[,i], at=order[i], ticks=median[i], twd=7, tlen=2, width=1/3)
}
denstrip(out2$sims.list$ICC, at=1, ticks=out2$q50$ICC, twd=7, tlen=2, width=1/3, colmax='red')
axis(1)
axis(2, las=1, at = 1:9,
    labels=c(expression(paste('Population growth rate (', italic(r), ')')),
        expression(paste('Double brooding (', kappa, ')')),
        expression(paste('Immigration (', omega, ')')),
        expression(paste('Fledging success 2 (', zeta[2], ')')),
        expression(paste('Fledging success 1 (', zeta[1], ')')),
        expression(paste('Clutch size 2 (', rho[2], ')')),
        expression(paste('Clutch size 1 (', rho[1], ')')),
        expression(paste('Adult survival (', phi[a], ')')),
        expression(paste('Juvenile survival (', phi[j], ')'))))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
