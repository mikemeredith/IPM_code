# Schaub & Kéry (2022) Integrated Population Models
# Chapter 7 : Assessment of integrated population models
# ------------------------------------------------------
# Code from final MS.

# Run time approx. 4 mins

library(IPMbook) ; library(jagsUI)

# ~~~ need these functions defined earlier ~~~
plotGOF <- function(jagsout, obs, rep, main=NA, showP=TRUE,
    ylab="Discrepancy replicate data", xlab="Discrepancy observed data",
    pch=16, cex = 0.8, col=1){
  OBS <- jagsout$sims.list[[obs]]
  REP <- jagsout$sims.list[[rep]]
  lim <- quantile(c(OBS, REP), c(0.0001, 0.999))
  plot(OBS, REP, pch=pch, cex=cex, ylim=lim, xlim=lim,
      ylab=ylab, xlab=xlab, main=main, axes=FALSE, col=col)
  axis(1); axis(2)
  segments(lim[1], lim[1], lim[2], lim[2], lty=3)
  bp <- round(mean(REP > OBS),2)
  if(showP){
    loc <- ifelse(bp < 0.5, "topleft", "bottomright")
    legend(loc, legend=bquote(p[B]==.(bp)), bty="n")
  }
  return(invisible(bp))
}

plotDS <- function(jagsname, param, param2, at, twd=c(3,1.5,1.5),
    tlen=c(2,2,2), width=1/3, colmax="black"){
  if(missing(param2))
    param2 <- param
  ticks <- eval(parse(text=paste0(jagsname, "$summary['", param2, "',c(1,3,7)]")))
  denstrip::denstrip(x=eval(parse(text=paste0(jagsname, "$sims.list$", param))),
      at=at, horiz=FALSE, ticks=ticks, twd=twd, tlen=tlen, width=width, colmax=colmax)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 7.2 Assumptions of integrated population models
# ===============================================

# 7.2.1 Assumptions made for the component data likelihoods
# ---------------------------------------------------------

# 7.2.1.4 Posterior predictive checks for IPMs with a hidden parameter
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# ~~~~ Models and fitting of the IPMs with fecundity as hidden parameter ~~~~
# 1. Small sample size

# Load data
library(IPMbook); library(jagsUI)
data(woodchat5)
marr <- marrayAge(woodchat5$ch, woodchat5$age)

# Bundle data
jags.data <- with(woodchat5, list(marr.j=marr[,,1], marr.a=marr[,,2],
    n.occasions=dim(marr)[2], rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
    C=count, pNinit=dUnif(1, 300)))

# IPM4: correctly specified
# Write JAGS model file
cat(file="model5.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  sigma ~ dunif(0.5, 1000)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois((N[1,t] + N[2,t]) * mean.f / 2 * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Assessing the fit of the state-space model
  # 1. Compute fit statistic for observed data
  # Discrepancy measure: mean absolute error
  for (t in 1:n.occasions){
    C.exp[t] <- N[1,t] + N[2,t]                            # Expected counts
    Dssm.obs[t] <- abs((C[t] - C.exp[t]) / C[t])           # Discrepancy measure
  }
  Dmape.obs <- sum(Dssm.obs)

  # 2. Compute fit statistic for replicate data
  # Discrepancy measure: mean absolute error
  for (t in 1:n.occasions){
    C.rep[t] ~ dnorm(N[1,t] + N[2,t], tau)                 # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t])   # Discrepancy measure
  }
  Dmape.rep <- sum(Dssm.rep)

  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1 - p[t]                                       # Probability of non-recapture
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

  # Assessing the fit of the capture-recapture model
  # 1. Compute fit statistic for observed data
  # Discrepancy measure: Freeman-Tukey statistic
  for (t in 1:(n.occasions-1)){
    for (j in 1:n.occasions){
      marr.j.exp[t,j] <- pr.j[t,j] * rel.j[t]              # Expected values
      marr.a.exp[t,j] <- pr.a[t,j] * rel.a[t]              # Expected values
      Dcjs.obs[t,j] <- pow(pow(marr.j[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.obs[t+n.occasions-1,j] <- pow(pow(marr.a[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.obs <- sum(Dcjs.obs)

  # 2. Compute fit statistic for replicate data
  # Discrepancy measure: Freeman-Tukey statistic
  for (t in 1:(n.occasions-1)){
    marr.j.rep[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])   # Generate replicate data
    marr.a.rep[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    for (j in 1:n.occasions){
      Dcjs.rep[t,j] <- pow(pow(marr.j.rep[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.rep[t+n.occasions-1,j] <- pow(pow(marr.a.rep[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.rep <- sum(Dcjs.rep)
}
")

# Initial values
inits <- function(){list(mean.p=runif(1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma",
    "Dmape.obs", "Dmape.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000

# Call JAGS (ART 1 min), check convergence and summarize posteriors
out9 <- jags(jags.data, inits, parameters, "model5.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out9)
print(out9, 3)



# IPM5: no age-dependence in survival
# Write JAGS model file
cat(file="model6.txt", "
model {
  # Priors and linear models
  mean.s ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.s
    sa[t] <- mean.s
    p[t] <- mean.p
  }


  sigma ~ dunif(0.5, 1000)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois((N[1,t] + N[2,t]) * mean.f / 2 * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Assessing the fit of the state-space model
  # 1. Compute fit statistic for observed data
  # Discrepancy measure: mean absolute error
  for (t in 1:n.occasions){
    C.exp[t] <- N[1,t] + N[2,t]                            # Expected counts
    Dssm.obs[t] <- abs((C[t] - C.exp[t]) / C[t])           # Discrepancy measure
  }
  Dmape.obs <- sum(Dssm.obs)

  # 2. Compute fit statistic for replicate data
  # Discrepancy measure: mean absolute error
  for (t in 1:n.occasions){
    C.rep[t] ~ dnorm(N[1,t] + N[2,t], tau)                 # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t])   # Discrepancy measure
  }
  Dmape.rep <- sum(Dssm.rep)

  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1 - p[t]                                       # Probability of non-recapture
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

  # Assessing the fit of the capture-recapture model
  # 1. Compute fit statistic for observed data
  # 1.1. Discrepancy measure: Freeman-Tukey statistic
  for (t in 1:(n.occasions-1)){
    for (j in 1:n.occasions){
      marr.j.exp[t,j] <- pr.j[t,j] * rel.j[t]              # Expected values
      marr.a.exp[t,j] <- pr.a[t,j] * rel.a[t]              # Expected values
      Dcjs.obs[t,j] <- pow(pow(marr.j[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.obs[t+n.occasions-1,j] <- pow(pow(marr.a[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.obs <- sum(Dcjs.obs)

  # 2. Compute fit statistic for replicate data
  # 2.1. Discrepancy measure: Freeman-Tukey statistic
  for (t in 1:(n.occasions-1)){
    marr.j.rep[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])   # Generate replicate data
    marr.a.rep[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
      for (j in 1:n.occasions){
      Dcjs.rep[t,j] <- pow(pow(marr.j.rep[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.rep[t+n.occasions-1,j] <- pow(pow(marr.a.rep[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.rep <- sum(Dcjs.rep)
  }
  ")

# Initial values
inits <- function(){list(mean.p=runif(1))}

# Parameters monitored
parameters <- c("mean.s", "mean.p", "mean.f", "N", "sigma",
    "Dmape.obs", "Dmape.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000

# Call JAGS (ART 1 min), check convergence and summarize posteriors
out10 <- jags(jags.data, inits, parameters, "model6.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out10)
print(out10, 3)

# IPM6: no reproduction of 1y shrikes
# Write JAGS model file
cat(file="model7.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  sigma ~ dunif(0.5, 1000)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois(N[2,t] * mean.f / 2 * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Assessing the fit of the state-space model
  # 1. Compute fit statistic for observed data
  # Discrepancy measure: mean absolute error
  for (t in 1:n.occasions){
    C.exp[t] <- N[1,t] + N[2,t]                            # Expected counts
    Dssm.obs[t] <- abs((C[t] - C.exp[t]) / C[t])           # Discrepancy measure
  }
  Dmape.obs <- sum(Dssm.obs)

  # 2. Compute fit statistic for replicate data
  # Discrepancy measure: mean absolute error
  for (t in 1:n.occasions){
    C.rep[t] ~ dnorm(N[1,t] + N[2,t], tau)                 # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t])   # Discrepancy measure
  }
  Dmape.rep <- sum(Dssm.rep)

  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1 - p[t]                                       # Probability of non-recapture
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

  # Assessing the fit of the capture-recapture model
  # 1. Compute fit statistic for observed data
  # 1.1. Discrepancy measure: Freeman-Tukey statistic
  for (t in 1:(n.occasions-1)){
    for (j in 1:n.occasions){
      marr.j.exp[t,j] <- pr.j[t,j] * rel.j[t]              # Expected values
      marr.a.exp[t,j] <- pr.a[t,j] * rel.a[t]              # Expected values
      Dcjs.obs[t,j] <- pow(pow(marr.j[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.obs[t+n.occasions-1,j] <- pow(pow(marr.a[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.obs <- sum(Dcjs.obs)

  # 2. Compute fit statistic for replicate data
  # 2.1. Discrepancy measure: Freeman-Tukey statistic
  for (t in 1:(n.occasions-1)){
    marr.j.rep[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])   # Generate replicate data
    marr.a.rep[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
      for (j in 1:n.occasions){
      Dcjs.rep[t,j] <- pow(pow(marr.j.rep[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.rep[t+n.occasions-1,j] <- pow(pow(marr.a.rep[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.rep <- sum(Dcjs.rep)
}
")

# Initial values
inits <- function(){list(mean.p=runif(1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma",
    "Dmape.obs", "Dmape.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000

# Call JAGS (ART 1 min), check convergence and summarize posteriors
out11 <- jags(jags.data, inits, parameters, "model7.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out11)
print(out11, 3)

# 2. Large sample size
# Analysis of the five times larger data set
# Load data
library(IPMbook)
data(woodchat7)
marr <- marrayAge(woodchat7$ch, woodchat7$age)

# Bundle data
jags.data <- with(woodchat7, list(marr.j=marr[,,1], marr.a=marr[,,2],
    n.occasions=dim(marr)[2], rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
    C=count, pNinit=dUnif(100, 1000)))

# Fit IPM4
# Initial values
inits <- function(){list(mean.p=runif(1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma",
    "Dmape.obs", "Dmape.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000

# Call JAGS (ART 8 min), check convergence and summarize posteriors
out12 <- jags(jags.data, inits, parameters, "model5.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out12)
print(out12, 3)

# Fit IPM5
# Initial values
inits <- function(){list(mean.p=runif(1))}

# Parameters monitored
parameters <- c("mean.s", "mean.p", "mean.f", "N", "sigma",
    "Dmape.obs", "Dmape.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000

# Call JAGS (ART 8 min), check convergence and summarize posteriors
out13 <- jags(jags.data, inits, parameters, "model6.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out13)
print(out13, 3)

# Fit IPM6
# Initial values
inits <- function(){list(mean.p=runif(1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma",
    "Dmape.obs", "Dmape.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000

# Call JAGS (ART 8 min), check convergence and summarize posteriors
out14 <- jags(jags.data, inits, parameters, "model7.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out14)
print(out14, 3)

save(out9, out10, out11, out12, out13, out14,
    file="ResultsChapter7.2.1.4.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for Fig. 7.6 ~~~~

load("ResultsChapter7.2.1.4.Rdata")

library(scales)
library(plotrix)
co <- viridis_pal(option='E')(20)[c(5, 11, 16)]

op <- par(mfrow=c(3,2), las=1, mar=c(3,7,3,1))
plotGOF(out9, "Dmape.obs", "Dmape.rep", main="State-space model", col=alpha(co[1], 0.3))
mtext(expression(bold(IPM[4])), side=2, las=0, line=5, font=2)
corner.label('A', font=2, cex=1.25)

par(mar=c(3,4,3,1))
plotGOF(out9, "DFT.obs", "DFT.rep", main="Cormack-Jolly-Seber model", ylab=NA, col=alpha(co[1], 0.3))
corner.label('B', font=2, cex=1.25)

par(mar=c(3,7,3,1))
plotGOF(out10, "Dmape.obs", "Dmape.rep", col=alpha(co[2], 0.3))
mtext(expression(bold(IPM[5])), side=2, las=0, line=5, font=2)
corner.label('C', font=2, cex=1.25)

par(mar=c(3,4,3,1))
plotGOF(out10, "DFT.obs", "DFT.rep", ylab=NA, col=alpha(co[2], 0.3))
corner.label('D', font=2, cex=1.25)

par(mar=c(5,7,3,1))
plotGOF(out11, "Dmape.obs", "Dmape.rep", col=alpha(co[3], 0.3))
mtext(expression(bold(IPM[6])), side=2, las=0, line=5, font=2)
corner.label('E', font=2, cex=1.25)

par(mar=c(5,4,3,1))
plotGOF(out11, "DFT.obs", "DFT.rep", ylab=NA, col=alpha(co[3], 0.3))
corner.label('F', font=2, cex=1.25)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for Fig. 7.7 ~~~~

load("ResultsChapter7.2.1.4.Rdata")

library(scales)
co <- viridis_pal(option='E')(20)[c(5, 11, 16)]

op <- par(mfrow=c(3,2), las=1, mar=c(3,7,3,1))
plotGOF(out12, "Dmape.obs", "Dmape.rep", main="State-space model", col=alpha(co[1], 0.3))
mtext(expression(bold(IPM[4])), side=2, las=0, line=5, font=2)
corner.label('A', font=2, cex=1.25)

par(mar=c(3,4,3,1))
plotGOF(out12, "DFT.obs", "DFT.rep", main="Cormack-Jolly-Seber model", ylab=NA, col=alpha(co[1], 0.3))
corner.label('B', font=2, cex=1.25)

par(mar=c(3,7,3,1))
plotGOF(out13, "Dmape.obs", "Dmape.rep", col=alpha(co[2], 0.3))
mtext(expression(bold(IPM[5])), side=2, las=0, line=5, font=2)
corner.label('C', font=2, cex=1.25)

par(mar=c(3,4,3,1))
plotGOF(out13, "DFT.obs", "DFT.rep", ylab=NA, col=alpha(co[2], 0.3))
corner.label('D', font=2, cex=1.25)

par(mar=c(5,7,3,1))
plotGOF(out14, "Dmape.obs", "Dmape.rep", col=alpha(co[3], 0.3))
mtext(expression(bold(IPM[6])), side=2, las=0, line=5, font=2)
corner.label('E', font=2, cex=1.25)

par(mar=c(5,4,3,1))
plotGOF(out14, "DFT.obs", "DFT.rep", ylab=NA, col=alpha(co[3], 0.3))
corner.label('F', font=2, cex=1.25)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for Fig. 7.8 ~~~~
load("ResultsChapter7.2.1.4.Rdata")

data(woodchat5)
data(woodchat7)

library(scales)
co <- viridis_pal(option='E')(20)[c(5, 11, 16)]

# Produce plot with counts against estimated population size from the 3 models
Ntot <- NtotL <- array(NA, dim=c(dim(out9$sims.list$N)[1], 3, dim(out9$sims.list$N)[3]))
for (t in 1:length(woodchat5$count)){
  Ntot[,1,t] <- out9$sims.list$N[,1,t] + out9$sims.list$N[,2,t]
  Ntot[,2,t] <- out10$sims.list$N[,1,t] + out10$sims.list$N[,2,t]
  Ntot[,3,t] <- out11$sims.list$N[,1,t] + out11$sims.list$N[,2,t]

  NtotL[,1,t] <- out12$sims.list$N[,1,t] + out12$sims.list$N[,2,t]
  NtotL[,2,t] <- out13$sims.list$N[,1,t] + out13$sims.list$N[,2,t]
  NtotL[,3,t] <- out14$sims.list$N[,1,t] + out14$sims.list$N[,2,t]
}
qu <- function(x) quantile(x, c(0.025, 0.975))
T <- length(woodchat5$count)
d <- 0.2

op <- par(mfrow=c(2,1), las=1, mar=c(3,5,3,1))
plot(y=apply(Ntot[,1,], 2, mean), x=(1:T)-d, type="b", ylim=c(80, 180), ylab="Number",
    xlab=NA, las=1, pch=16, axes=FALSE, col=co[1], main="Small sample size")
axis(2)
axis(1, at=1:20, tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20), tcl=-0.5, labels=NA)

segments((1:T)-d, apply(Ntot[,1,], 2, qu)[1,], (1:T)-d, apply(Ntot[,1,], 2, qu)[2,], col=co[1])
points(y=apply(Ntot[,2,], 2, mean), x=1:T, type="b", pch=16, col=co[2])
segments(1:T, apply(Ntot[,2,], 2, qu)[1,], 1:T, apply(Ntot[,2,], 2, qu)[2,], col=co[2])
points(y=apply(Ntot[,3,], 2, mean), x=(1:T)+d, type="b", pch=16, col=co[3])
segments((1:T)+d, apply(Ntot[,3,], 2, qu)[1,], (1:T)+d, apply(Ntot[,3,], 2, qu)[2,], col=co[3])
points(woodchat5$count, type="b", pch=1, lty=3)

par(mar=c(5,5,1,1))
plot(y=apply(NtotL[,1,], 2, mean), x=(1:T)-d, type="b", ylim=c(400, 900), ylab="Number",
    xlab="Year", las=1, pch=16, axes=FALSE, col=co[1], main="Large sample size")
axis(2)
axis(1, at=1:20, tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20), tcl=-0.5, labels=c(5, 10, 15, 20))

segments((1:T)-d, apply(NtotL[,1,], 2, qu)[1,], (1:T)-d, apply(NtotL[,1,], 2, qu)[2,], col=co[1])
points(y=apply(NtotL[,2,], 2, mean), x=1:T, type="b", pch=16, col=co[2])
segments(1:T, apply(NtotL[,2,], 2, qu)[1,], 1:T, apply(NtotL[,2,], 2, qu)[2,], col=co[2])
points(y=apply(NtotL[,3,], 2, mean), x=(1:T)+d, type="b", pch=16, col=co[3])
segments((1:T)+d, apply(NtotL[,3,], 2, qu)[1,], (1:T)+d, apply(NtotL[,3,], 2, qu)[2,], col=co[3])
points(woodchat7$count, type="b", pch=1, lty=3)
legend("topleft", legend=c("Observed counts", expression("Estimates from IPM"[4]),
    expression("Estimates from IPM"[5]), expression("Estimates from IPM"[6])),
    pch=c(1, rep(16, 3)), col=c("black", co), lty=c(3, 1, 1, 1), bty="n")
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for Fig. 7.9 ~~~~

load("ResultsChapter7.2.1.4.Rdata")

library(denstrip)
library(scales)
co <- viridis_pal(option='E')(20)[c(5, 11, 16)]
op <- par(las=1)
plot(0, ylim=c(2, 7.3), xlim=c(0.7,7.3), axes=FALSE, pch=NA, xlab=NA, ylab="Productivity")
plotDS("out9", "mean.f", at=1, colmax=co[1])
plotDS("out10", "mean.f", at=2, colmax=co[2])
plotDS("out11", "mean.f", at=3, colmax=co[3])
plotDS("out12", "mean.f", at=5, colmax=co[1])
plotDS("out13", "mean.f", at=6, colmax=co[2])
plotDS("out14", "mean.f", at=7, colmax=co[3])
axis(2)
labs <- c(expression(IPM[4]), expression(IPM[5]), expression(IPM[6]),
    expression(IPM[4]), expression(IPM[5]), expression(IPM[6]))
axis(1, at=c(1:3, 5:7), labels=labs)
mtext("Small sample size", side=1, at=2, line=3)
mtext("Large sample size", side=1, at=6, line=3)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
