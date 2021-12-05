# Schaub & Kéry (2022) Integrated Population Models
# Chapter 5 : Introduction to integrated population models
# --------------------------------------------------------

# Run time approx. 2 mins

library(IPMbook) ; library(jagsUI)

# ~~~ requires data prepared in section 5.2 ~~~
data(woodchat5)
marr <- marrayAge(woodchat5$ch, woodchat5$age)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5.3 Our first integrated population model
# =========================================

# Bundle data and produce data overview
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=woodchat5$repro[,1],
    age=woodchat5$repro[,3], C=woodchat5$count)
str(jags.data)
# List of 8
# $ marr.j     : num [1:19, 1:20] 8 0 0 0 0 0 0 0 0 0 ...
# $ marr.a     : num [1:19, 1:20] 16 0 0 0 0 0 0 0 0 0 ...
# $ n.occasions: int 20
# $ rel.j      : num [1:19] 51 53 55 65 73 66 61 76 65 75 ...
# $ rel.a      : num [1:19] 36 39 44 61 61 50 43 61 51 53 ...
# $ J          : num [1:929] 6 2 2 5 3 5 3 2 3 2 ...
# $ age        : num [1:929] 1 1 1 1 1 1 1 1 1 1 ...
# $ C          : num [1:20] 91 119 131 88 139 145 148 116 112 106 ...

# Write JAGS model file
cat(file="model3.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f[1] ~ dunif(0, 10)
  mean.f[2] ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial stage-specific population sizes: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- N[1,t] * mean.f[1] / 2 * mean.sj + N[2,t] * mean.f[2] / 2 * mean.sj
    N[2,t+1] <- (N[1,t] + N[2,t]) * mean.sa
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Productivity data (Poisson regression model)
  for (i in 1:length(J)){
    J[i] ~ dpois(mean.f[age[i]])
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
    q[t] <- 1 - p[t]                    # Probability of non-recapture
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
  # Annual population growth rate
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])
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
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma", "ann.growth.rate", "Ntot")

# MCMC settings
ni <- 40000; nb <- 10000; nc <- 3; nt <- 3; na <- 3000

# Call JAGS (ART 2 min), check convergence and summarize posteriors
out3 <- jags(jags.data, inits, parameters, "model3.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    # n.thin=nt, n.adapt=na)
    n.thin=nt, n.adapt=na, parallel=TRUE)  # ~~~ for testing
traceplot(out3)
print(out3, 3)
#                         mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# mean.sj                0.298  0.010    0.279    0.298    0.317    FALSE 1 1.002  1175
# mean.sa                0.539  0.012    0.515    0.539    0.564    FALSE 1 1.003   857
# mean.p                 0.608  0.019    0.571    0.608    0.645    FALSE 1 1.000  9184
# mean.f[1]              2.671  0.076    2.526    2.670    2.822    FALSE 1 1.000  7172
# mean.f[2]              3.674  0.085    3.510    3.672    3.844    FALSE 1 1.001  4281
# N[1,1]                47.751 32.297    2.676   43.470  110.762    FALSE 1 1.006   371
# N[2,1]                62.412 29.082    6.701   65.476  106.317    FALSE 1 1.006   380
# N[1,2]                53.101  4.429   44.608   53.092   61.675    FALSE 1 1.004   499
# N[2,2]                59.402  4.247   51.433   59.259   68.046    FALSE 1 1.002   897
# [ ...output truncated... ]
# N[1,20]               70.264  4.264   61.993   70.211   78.927    FALSE 1 1.000 12067
# N[2,20]               79.456  4.399   71.016   79.381   88.332    FALSE 1 1.001  2774
# sigma                 16.976  3.149   12.116   16.516   24.414    FALSE 1 1.000 10870
# ann.growth.rate[1]     1.023  0.041    0.945    1.026    1.087    FALSE 1 1.006   363
# ann.growth.rate[2]     1.016  0.007    1.003    1.015    1.029    FALSE 1 1.002   979
# [ ...output truncated... ]
# ann.growth.rate[18]    1.016  0.005    1.006    1.016    1.027    FALSE 1 1.000 18640
# ann.growth.rate[19]    1.016  0.005    1.006    1.016    1.027    FALSE 1 1.000 18640
# Ntot[1]              110.163  7.608   95.745  109.918  125.721    FALSE 1 1.001  1648
# Ntot[2]              112.503  6.546   99.570  112.493  125.646    FALSE 1 1.000 11582
# [ ...output truncated... ]
# Ntot[19]             147.325  7.114  133.524  147.263  161.763    FALSE 1 1.000 30000
# Ntot[20]             149.720  7.854  134.538  149.609  165.765    FALSE 1 1.000 30000


# ~~~~ code for Figure 5.2 ~~~~
mag <- 1.25
cex.tif <- mag * 1.25
lwd.tif <- 3 * mag
op <- par(mar=c(4, 4, 3, 0), las=1, cex=cex.tif, lwd=lwd.tif)
u <- col2rgb("grey82")
T <- length(woodchat5$count)
col.pol <- rgb(u[1], u[2], u[3], alpha=100, maxColorValue=255)
plot(out3$mean$Ntot, type="n",
    ylim=range(c(out3$q2.5$Ntot, out3$q97.5$Ntot, woodchat5$count)),
    ylab="Number", xlab="Year", las=1, cex=1.5, axes=FALSE)
axis(2, las=1, lwd=lwd.tif)
axis(2, at=c(90, 110, 130, 150), labels=NA, tcl=-0.25, lwd=lwd.tif)
axis(1, at=1:T, labels=NA, tcl=-0.25, lwd=lwd.tif)
axis(1, at=c(5, 10, 15, 20), labels=c(5, 10, 15, 20), tcl=-0.5, lwd=lwd.tif)
polygon(c(1:T, T:1), c(out3$q2.5$Ntot, out3$q97.5$Ntot[T:1]), border=NA, col=col.pol)
points(out3$mean$Ntot, type="b", col="black", pch=16, lty=1, lwd=lwd.tif)
points(woodchat5$count, type="b", col="blue", pch=1, lty=2, lwd=lwd.tif)
legend("topleft", legend=c("Observed population counts", "Estimated population size"),
    pch=c(1, 16), lwd=c(lwd.tif, lwd.tif), col=c("blue", "black"), lty=c(2, 1), bty="n")
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
