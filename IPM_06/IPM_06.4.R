# Schaub & Kéry (2022) Integrated Population Models
# Chapter 6 : Benefits of integrated population modeling
# ------------------------------------------------------

# Run time approx. 5 mins

# 6.4 Estimation of process correlation among demographic parameters
# ==================================================================

library(IPMbook); library(jagsUI)
data(woodchat64)
marr <- marrayAge(woodchat64$ch, woodchat64$age)

# Bundle data
jags.data<- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C=woodchat64$count, J=woodchat64$J,
    B=woodchat64$B)

# Write JAGS model file
cat(file="model4.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.occasions-1)){
    logit(sj[t]) <- eps[1,t]
    logit(sa[t]) <- eps[2,t]
    log(f[t]) <- eps[3,t]
    eps[1:3,t] ~ dmnorm.vcov(mu[], sigma2[,])
    p[t] <- mean.p
  }

  mu[1] <- logit(mean.sj)
  mu[2] <- logit(mean.sa)
  mu[3] <- log(mean.f)

  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)
  mean.p ~ dunif(0, 1)

  # Reparameterize sigma2 in terms of SD and rho
  for (i in 1:3){
    sigma2[i,i] <- sigma[i] * sigma[i]
  }
  for (i in 1:2){
    for (j in (i+1):3){
      sigma2[i,j] <- sigma[i] * sigma[j] * rho[i+j-2]
      sigma2[j,i] <- sigma2[i,j]
    } #j
  } #i

  # Specify priors for SD and rho
  for (i in 1:3){
    sigma[i] ~ dunif(0, 5)
    rho[i] ~ dunif(-1,1)
  }

  sigma.res ~ dunif(0.5, 100)
  tau.res <- pow(sigma.res, -2)

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time: our model for population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- f[t] / 2 * sj[t] * (N[1,t] + N[2,t])
    N[2,t+1] <- sa[t] * (N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau.res)
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
    q[t] <- 1 - p[t] # Probability of non-recapture
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
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "sj", "sa", "f", "N", "sigma.res",
    "sigma2", "rho")

# MCMC settings
ni <- 60000; nb <- 10000; nc <- 3; nt <- 10; na <- 2000

# Call JAGS (ART 11 min), check convergence and summarize posteriors
out4 <- jags(jags.data, inits, parameters, "model4.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out4)
print(out4, 3)
#                mean     sd    2.5%     50%   97.5% overlap0     f  Rhat n.eff
# [ ...output truncated... ]
# sigma.res    11.305  4.138   3.261  11.057  20.272    FALSE 1.000 1.004   563
# sigma2[1,1]   0.044  0.031   0.005   0.037   0.123    FALSE 1.000 1.003   684
# sigma2[2,1]   0.019  0.017  -0.009   0.016   0.059     TRUE 0.912 1.002  1551
# sigma2[3,1]   0.001  0.015  -0.030   0.001   0.030     TRUE 0.534 1.003   839
# sigma2[1,2]   0.019  0.017  -0.009   0.016   0.059     TRUE 0.912 1.002  1551
# sigma2[2,2]   0.045  0.036   0.004   0.037   0.139    FALSE 1.000 1.008   395
# sigma2[3,2]   0.017  0.017  -0.011   0.015   0.054     TRUE 0.890 1.004   832
# sigma2[1,3]   0.001  0.015  -0.030   0.001   0.030     TRUE 0.534 1.003   839
# sigma2[2,3]   0.017  0.017  -0.011   0.015   0.054     TRUE 0.890 1.004   832
# sigma2[3,3]   0.047  0.022   0.020   0.043   0.101    FALSE 1.000 1.001 15000
# rho[1]        0.486  0.324  -0.250   0.535   0.948     TRUE 0.912 1.005   420
# rho[2]        0.030  0.307  -0.548   0.029   0.637     TRUE 0.534 1.004   582
# rho[3]        0.406  0.309  -0.267   0.440   0.905     TRUE 0.890 1.008   255

round(out4$mean$sigma2, 3)              # Posterior means of var-covar matrix
#       [,1]  [,2]  [,3]
# [1,] 0.044 0.019 0.001
# [2,] 0.019 0.045 0.017
# [3,] 0.001 0.017 0.047

# ~~~~ code for Figure 6.6 ~~~~
lwd.tif <- 2
plot(density(out4$sims.list$rho[,1]), xlim=c(-1, 1), axes=FALSE, main=NA,
    ylab="Density", xlab="Process correlation", lwd=lwd.tif)
lines(density(out4$sims.list$rho[,2]), col="red", lwd=lwd.tif)
lines(density(out4$sims.list$rho[,3]), col="blue", lwd=lwd.tif)
axis(1)
axis(2, las=1)
legend("topleft", lwd = c(lwd.tif, 3), col = c("red", "blue", "black"), bty = "n",
    legend = c(expression(rho[italic(s)[italic(j)]][italic(f)]),
    expression(rho[italic(s)[italic(a)]][italic(f)]),
    expression(rho[italic(s)[italic(j)]][italic(s)[italic(a)]])))
abline(v = 0, lwd = 1, lty = 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
