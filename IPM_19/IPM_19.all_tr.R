# Schaub & Kéry (2022) Integrated Population Models
# Chapter 19 : Grey catbird
# -------------------------

# Run time for test script 17 mins, full run 3.5 hrs

# 19.4 Component data likelihoods
# ===============================

library(IPMbook); library(jagsUI)
data(catbird)
str(catbird)
# List of 9
# $ y       : num [1:4276, 1:17] 0 1 0 0 1 1 0 0 0 0 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:4276] "427696" "427697" "427699" "427703" ...
# .. ..$ : chr [1:17] "1992" "1993" "1994" "1995" ...
# $ r       : num [1:4276] 0 0 0 0 0 0 0 0 1 0 ...
# $ station : num [1:4276] 1 1 1 1 1 1 1 1 1 1 ...
# $ count   : num [1:1298] 23 10 14 39 22 26 28 30 23 27 ...
# $ stratum : num [1:1298] 1 1 1 1 1 1 1 1 1 1 ...
# $ year    : num [1:1298] 1992 1992 1992 1992 1992 ...
# $ observer: num [1:1298] 1 2 3 4 5 6 7 8 9 8 ...
# $ firstyr : num [1:1298] 0 1 0 0 0 0 0 0 1 0 ...
# $ area    : num [1:9] 10982 14868 14329 2132 4192 ...


# 19.4.1 Population count data (BBS data) (no code)
# 19.4.2 Capture-recapture data (MAPS data) (no code)

# 19.5 The integrated population model
# ====================================

# Bundle data and produce data overview
jags.data <- with(catbird, list(y=y, f=getFirst(y), r=r, sta=station,
    n.sta=length(unique(station)), n.years=ncol(y), n.ind=nrow(y), C=count, n.counts=length(count),
    str=stratum, yr=year-1991, n.str=length(unique(stratum)), I=firstyr, obs=observer,
    n.obs=length(unique(observer)), area=area, pNinit=dUnif(1, 50)))
str(jags.data)
# List of 17
# $ y       : num [1:4276, 1:17] 0 1 0 0 1 1 0 0 0 0 ...
# $ f       : Named int [1:4276] 4 1 4 4 1 1 4 4 4 4 ...
# $ r       : num [1:4276] 0 0 0 0 0 0 0 0 1 0 ...
# $ sta     : num [1:4276] 1 1 1 1 1 1 1 1 1 1 ...
# $ n.sta   : int 38
# $ n.years : int 17
# $ n.ind   : int 4276
# $ C       : num [1:1298] 23 10 14 39 22 26 28 30 23 27 ...
# $ n.counts: int 1298
# $ str     : num [1:1298] 1 1 1 1 1 1 1 1 1 1 ...
# $ yr      : num [1:1298] 1 1 1 1 1 1 1 1 1 1 ...
# $ n.str   : int 9
# $ I       : num [1:1298] 0 1 0 0 0 0 0 0 1 0 ...
# $ obs     : num [1:1298] 1 2 3 4 5 6 7 8 9 8 ...
# $ n.obs   : int 172
# $ area    : num [1:9] 10982 14868 14329 2132 4192 ...
# $ pNinit  : num [1:50] 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 ...

# Write JAGS model file
cat(file="model1.txt", "
model {
  # Priors and linear models
  for (i in 1:n.ind){
    for (t in f[i]:(n.years-1)){
      logit(p[i,t]) <- lp0[t] + alpha[sta[i]]
    } #t
  } #i

  for (t in 1:(n.years-1)){
    phi[t] ~ dunif(0, 1)                                  # Survival probability
    lp0[t] <- logit(p0[t])                                # Recapture probability
    p0[t] ~ dunif(0, 1)
    kappa[t] ~ dunif(0, 1)                                # Predetermined residency probability
    pi[t] ~ dunif(0, 1)                                   # Residency probability
    rho[t] ~ dunif(0, 5)                                  # Recruitment
  }

  # Random station effect for recapture probability
  for (j in 1:n.sta){
    alpha[j] ~ dnorm(0, tau.p)
  }
  sigma.p ~ dunif(0, 10)
  tau.p <- pow(sigma.p, -2)

  # Random observer effect for survey observation model
  for (i in 1:n.obs){
    omega[i] ~ dnorm(0, tau.om)
  }
  tau.om ~ dgamma(0.001,0.001)
  sd.om <- pow(tau.om, -0.5)

  # Overdispersion
  for (k in 1:n.counts){
    eps[k] ~ dnorm(0, tau.eps)
  }
  tau.eps ~ dgamma(0.001,0.001)
  sd.eps <- pow(tau.eps, -0.5)

  # Start-up (novice) effect for survey
  eta ~ dnorm(0, 0.001)

  # BBS data (state-space model)
  # Model for initial abundance in each stratum: uniform priors
  for (s in 1:n.str){
    N[s,1] ~ dcat(pNinit)
  }
  # Process model over space and time: our model of population dynamics
  for (s in 1:n.str){
    for (t in 2:n.years){
      S[s,t] ~ dbin(phi[t-1], N[s,t-1])
      R[s,t] ~ dpois(rho[t-1] * N[s,t-1])
      N[s,t] <- S[s,t] + R[s,t]
    } #t
  } #s
  # Observation model: overdispersed Poisson
  for (k in 1:n.counts){
    log(lambda[k]) <- log(N[str[k],yr[k]]) + omega[obs[k]] + eta * I[k] + eps[k]
    C[k] ~ dpois(lambda[k])
  }

  # MAPS data (CJS model with state-space likelihood)
  for (i in 1:n.ind){
    # Model for initial state
    for (t in 1:f[i]){
      z[i,t] <- 1
    } #t
    # Model for residency state
    u[i] ~ dbern(pi[f[i]])
    r[i] ~ dbern(u[i] * kappa[f[i]])                      # Observation model for residency
    # State process
    for (t in (f[i]+1):n.years){
      z[i,t] ~ dbern(z[i,t-1] * phi[t-1] * u[i])
      # Observation process
      y[i,t] ~ dbern(p[i,t-1] * z[i,t])
    } #t
  } #i

  # Derived parameters
  for (s in 1:n.str){
    for (t in 1:n.years){
      wN[s,t] <- area[s] * N[s,t] / sum(area)
    } #t
    trend[s] <- pow(N[s,n.years] / N[s,1], 1/(n.years-1)) # Strata-spe.trend
  } #s
  for (t in 1:n.years){
    wIndex[t] <- sum(wN[1:n.str,t])                       # Weighted abundance index
  }
  trend.com <- pow(wIndex[n.years] / wIndex[1], 1/(n.years-1)) # Combined trend
}
")

# Initial values
ny <- jags.data$n.years
ns <- jags.data$n.str
uini <- rep(1, nrow(catbird$y))
inits <- function(){
  Nini <- matrix(NA, nrow=ns, ncol=ny)
  Nini[,1] <- rpois(ns, 17)
  list(z=zInit(catbird$y), u=uini, N=Nini, rho=runif(ny-1, 0.3, 1.3))
}

# Parameters monitored
parameters <- c("phi", "rho", "pi", "p0", "kappa", "eta", "alpha", "sd.om", "sd.eps", "wIndex", "trend",
    "trend.com", "N")

# MCMC settings
# ni <- 70000; nb <- 20000; nc <- 3; nt <- 25; na <- 2000
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 2000  # ~~~ for testing, 13 mins

# Call JAGS from R (ART 191 min) and check convergence
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1)

# ~~~ save the results ~~~
save(out1, file="CatbirdResults.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~

# 19.6 Results
# ============

print(out1, 3)
                   # mean      sd      2.5%       50%     97.5% overlap0     f  Rhat n.eff
# phi[1]            0.804   0.137     0.484     0.828     0.991    FALSE 1.000 1.001  2576
# phi[2]            0.719   0.157     0.409     0.727     0.978    FALSE 1.000 1.002   987
# [... output truncated ...]
# phi[15]           0.561   0.151     0.323     0.539     0.916    FALSE 1.000 1.004   511
# phi[16]           0.707   0.211     0.269     0.759     0.985    FALSE 1.000 1.017   153
# rho[1]            0.266   0.161     0.029     0.245     0.629    FALSE 1.000 1.001  2296
# rho[2]            0.280   0.177     0.017     0.265     0.659    FALSE 1.000 1.003   916
# [... output truncated ...]
# rho[15]           0.484   0.187     0.098     0.493     0.830    FALSE 1.000 1.004   594
# rho[16]           0.284   0.223     0.010     0.224     0.777    FALSE 1.000 1.016   151
# pi[1]             0.322   0.075     0.203     0.312     0.496    FALSE 1.000 1.000  6000
# pi[2]             0.519   0.106     0.354     0.507     0.770    FALSE 1.000 1.004   569
# [... output truncated ...]
# pi[15]            0.405   0.088     0.262     0.396     0.604    FALSE 1.000 1.001  5991
# pi[16]            0.538   0.126     0.335     0.524     0.825    FALSE 1.000 1.001  2668
# p0[1]             0.494   0.118     0.269     0.492     0.728    FALSE 1.000 1.002  1002
# p0[2]             0.301   0.083     0.163     0.294     0.483    FALSE 1.000 1.000  6000
# [... output truncated ...]
# p0[15]            0.309   0.090     0.158     0.300     0.505    FALSE 1.000 1.002  2213
# p0[16]            0.296   0.163     0.117     0.245     0.764    FALSE 1.000 1.018   171
# kappa[1]          0.372   0.098     0.198     0.366     0.573    FALSE 1.000 1.000  6000
# kappa[2]          0.304   0.071     0.179     0.299     0.450    FALSE 1.000 1.003   755
# [... output truncated ...]
# kappa[15]         0.236   0.061     0.133     0.230     0.371    FALSE 1.000 1.000  6000
# kappa[16]         0.283   0.075     0.158     0.276     0.453    FALSE 1.000 1.002  1087
# eta               0.008   0.054    -0.096     0.008     0.112     TRUE 0.561 1.001  2773
# alpha[1]         -0.063   0.245    -0.529    -0.070     0.417     TRUE 0.607 1.000  3298
# alpha[2]         -0.111   0.731    -1.567    -0.103     1.325     TRUE 0.561 1.000  6000
# [... output truncated ...]
# alpha[37]        -0.199   0.713    -1.687    -0.177     1.192     TRUE 0.607 1.001  2609
# alpha[38]        -0.098   0.718    -1.552    -0.075     1.314     TRUE 0.550 1.001  6000
# sd.om             0.672   0.046     0.589     0.669     0.768    FALSE 1.000 1.000  6000
# sd.eps            0.294   0.014     0.268     0.294     0.321    FALSE 1.000 1.000  3926
# wIndex[1]        19.198   1.498    16.393    19.158    22.287    FALSE 1.000 1.003   850
# wIndex[2]        20.443   1.586    17.508    20.357    23.706    FALSE 1.000 1.003   747
# [... output truncated ...]
# wIndex[16]       20.775   1.694    17.676    20.712    24.341    FALSE 1.000 1.002  1261
# wIndex[17]       20.448   1.775    17.177    20.363    24.195    FALSE 1.000 1.003   795
# trend[1]          1.020   0.011     0.998     1.020     1.042    FALSE 1.000 1.000  6000
# trend[2]          0.998   0.009     0.981     0.998     1.016    FALSE 1.000 1.000  1846
# [... output truncated ...]
# trend[8]           1.010   0.023    0.965     1.010     1.055    FALSE 1.000 1.000  4797
# trend[9]           0.961   0.028    0.917     0.958     1.000    FALSE 1.000 1.000  1245
# trend.com          1.004   0.006    0.9993    1.004     1.015    FALSE 1.000 1.000  6000
# N[1,1]           25.999   4.257    19.000    26.000    35.000    FALSE 1.000 1.003   717
# N[2,1]           23.695   3.433    18.000    23.000    31.000    FALSE 1.000 1.009   280
# [... output truncated ...]
# N[8,17]          36.925  14.075    16.000    35.000    70.000    FALSE 1.000 1.001  4363
# N[9,17]           1.082   0.291     1.000     1.000     2.000    FALSE 1.000 1.000  6000

# ~~~~ code for Fig. 19.4 ~~~~
op <- par(mfrow=c(3,1), mar=c(2.5, 4.5, 0.5, 1))

years <- 1992:2008
ny <- length(years)
ti <- seq(from=1.5, by=1, length.out=ny-1)

plot(x=ti, y=out1$mean$phi, type="b", pch=16, ylim=c(0, 1), ylab="Apparent survival",
    xlab=NA, axes=FALSE)
segments(ti, out1$q2.5$phi, ti, out1$q97.5$phi)
axis(1, at=1:ny, labels=NA, tcl=-0.25)
axis(1, at=c(1, 4, 7, 10, 13, 16), labels=NA, tcl=-0.5)
axis(2, las=1)

par(mar=c(2.5, 4.5, 0.5, 1))
plot(x=ti, y=out1$mean$rho, type="b", pch=16, ylim=c(0, 1.1), ylab="Recruitment",
    xlab=NA, axes=FALSE)
segments(ti, out1$q2.5$rho, ti, out1$q97.5$rho)
axis(1, at=1:ny, labels=NA, tcl=-0.25)
axis(1, at=c(1, 4, 7, 10, 13, 16), labels=NA, tcl=-0.5)
axis(2, las=1)

par(mar=c(4.5, 4.5, 0.5, 1))
plot(x=1:ny, y=out1$mean$wIndex, type="b", pch=16, ylim=c(15, 30), ylab="Abundance index",
    xlab=NA, axes=FALSE)
segments(1:ny, out1$q2.5$wIndex, 1:ny, out1$q97.5$wIndex)
axis(1, at=1:ny, labels=NA, tcl=-0.25)
axis(1, at=c(1, 4, 7, 10, 13, 16), labels=years[c(1, 4, 7, 10, 13, 16)], tcl=-0.5)
axis(2, las=1)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 19.5 ~~~~
library(denstrip)

op <- par(mar=c(3.5, 5, 1, 1))
plot(0, xlim=c(0.8, 10.2), ylim=c(0.9,1.1), axes=FALSE, pch=NA, ylab='Abundance trend', xlab=NA)
for (s in 1:9){
  denstrip(out1$sims.list$trend[,s], at=s,
      ticks=c(out1$mean$trend[s], out1$q2.5$trend[s], out1$q97.5$trend[s]),
      twd=c(5,2.5,2.5), tlen=c(2,2,2), width=1/5, horiz=FALSE)
   }
denstrip(out1$sims.list$trend.com, at=10,
    ticks=c(out1$mean$trend.com, out1$q2.5$trend.com, out1$q97.5$trend.com),
    twd=c(5,2.5,2.5), tlen=c(2,2,2), width=1/5, horiz=FALSE, colmax='red')
abline(h=1, lty=2)
axis(1, at=1:9, labels=1:9)
axis(2, las=1)
mtext('Stratum', side=1, line=2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 19.6 ~~~~
library(scales)
data(catbird)

years <- 1992:2008
ny <- length(years)

# Calculate statistics for the counts (across routes, per strata and year)
y <- with(catbird, cbind(count, stratum, year-1991))
count.stat <- array(NA, dim=c(9, 17, 7))
for (s in 1:9){
  for (t in 1:17){
    z <- which(y[,2] == s & y[,3] == t)
    count.stat[s,t,1] <- mean(y[z,1])
    count.stat[s,t,2] <- min(y[z,1])
    count.stat[s,t,3] <- max(y[z,1])
    count.stat[s,t,4] <- median(y[z,1])
    count.stat[s,t,5] <- quantile(y[z,1], 0.1)
    count.stat[s,t,6] <- quantile(y[z,1], 0.9)
    } #t
  } #s
count.stat[,,7] <- table(catbird$stratum, catbird$year)
count.stat[which(is.infinite(count.stat[,,]))] <- NA

op <- par(mfrow=c(3, 3), mar=c(2.5, 4.5, 1, 1))
for(s in 1:9){
  lims <- range(c(count.stat[s,,5], count.stat[s,,6], out1$q2.5$N[s,], out1$q97.5$N[s,]),
      na.rm=TRUE) * c(0.9, 1.1)
  plot(x=1:ny, y=out1$mean$N[s,], type="b", pch=16, ylim=lims, ylab="Abundance index",
      xlab=NA, axes=FALSE)
  if(s!=4 & s!=8)
    polygon(c(1:ny, ny:1), c(count.stat[s,,5], rev(count.stat[s,,6])),
        col=alpha('red', alpha=0.15), border=NA)
  if(s==4){
    polygon(c(1:7, 7:1), c(count.stat[s,1:7,5], rev(count.stat[s,1:7,6])),
        col=alpha('red', alpha=0.15), border=NA)
    polygon(c(9:ny, ny:9), c(count.stat[s,9:ny,5], rev(count.stat[s,9:ny,6])),
        col=alpha('red', alpha=0.15), border=NA)
  }
  if(s==8){
    polygon(c(1:5, 5:1), c(count.stat[s,1:5,5], rev(count.stat[s,1:5,6])),
        col=alpha('red', alpha=0.15), border=NA)
    polygon(c(7:ny, ny:7), c(count.stat[s,7:ny,5], rev(count.stat[s,7:ny,6])),
        col=alpha('red', alpha=0.15), border=NA)
  }

  segments(1:ny, out1$q2.5$N[s,], 1:ny, out1$q97.5$N[s,])
  lines(count.stat[s,,4], col="red")
  points(x=1:ny, y=out1$mean$N[s,], type="b", pch=16)
  axis(1, at=1:ny, labels=NA, tcl=-0.25)
  axis(1, at=c(1, 4, 7, 10, 13, 16), labels=years[c(1, 4, 7, 10, 13, 16)], tcl=-0.5)
  axis(2, las=1)
  text(y=lims[2], x=4.5, paste('stratum',s,': n = ', round(mean(count.stat[s,,7]),1)))
}
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
