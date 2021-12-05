# Schaub & Kéry (2022) Integrated Population Models
# Chapter 7 : Assessment of integrated population models
# ------------------------------------------------------

# Run time approx. 8 mins

# 7.2 Assumptions of integrated population models
# ===============================================

# 7.2.1 Assumptions made for the component data likelihoods
# ---------------------------------------------------------

# 7.2.1.1 Principles of posterior predictive checks (no code)

# 7.2.1.2 Application of posterior predictive checks
# ''''''''''''''''''''''''''''''''''''''''''''''''''

# Load data
library(IPMbook); library(jagsUI)
data(woodchat5)
marr <- marrayAge(woodchat5$ch, woodchat5$age)

# Bundle data
jags.data <- with(woodchat5, list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=repro[,1], year=repro[,2], age=repro[,3],
    C=count, pNinit=dUnif(1, 300)))
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
cat(file="model1.txt", "
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

  for (t in 1:n.occasions){
    f[1,t] <- mean.f[1]
    f[2,t] <- mean.f[2]
  }

  sigma ~ dunif(0.5, 1000)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
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

  # Assessing the fit of the state-space model
  # 1. Compute fit statistic for observed data
  # 1.1. Discrepancy measure: mean absolute error
  for (t in 1:n.occasions){
    C.exp[t] <- N[1,t] + N[2,t]                       # Expected counts
    Dssm.obs[t] <- abs((C[t] - C.exp[t]) / C[t])      # Discrepancy measure
  }
  Dmape.obs <- sum(Dssm.obs)
  # 1.2. Test statistic: number of turns
  for (t in 1:(n.occasions-2)){
    Tt1.obs[t] <- step(C[t+2] - C[t+1])
    Tt2.obs[t] <- step(C[t+1] - C[t])
    Tt3.obs[t] <- equals(Tt1.obs[t] + Tt2.obs[t], 1)
  }
  Tturn.obs <- sum(Tt3.obs)

  # 2. Compute fit statistic for replicate data
  # 2.1. Discrepancy measure: mean absolute error
  for (t in 1:n.occasions){
    C.rep[t] ~ dnorm(N[1,t] + N[2,t], tau)                # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t])  # Discrepancy measure
  }
  Dmape.rep <- sum(Dssm.rep)
  # 2.2. Test statistic: number of turns
  for (t in 1:(n.occasions-2)){
    Tt1.rep[t] <- step(C.rep[t+2] - C.rep[t+1])
    Tt2.rep[t] <- step(C.rep[t+1] - C.rep[t])
    Tt3.rep[t] <- equals(Tt1.rep[t] + Tt2.rep[t], 1)
  }
  Tturn.rep <- sum(Tt3.rep)

  # Productivity data (Poisson regression model)
  for (i in 1:length(J)){
    J[i] ~ dpois(f[age[i], year[i]])
  }

  # Assessing the fit of the Poisson regression model
  # 1. Compute fit statistic for observed data
  # 1.1. Discrepancy measure: deviance
  for (i in 1:length(J)){
    J.exp[i] <- f[age[i],year[i]]                     # Expected data
    D.obs[i] <- J[i] * log(J[i]/J.exp[i]) - (J[i] - J.exp[i])
  }
  Dd.obs <- sum(D.obs)
  # 1.2. Compute test statistic variance-mean ratio
  Tvm.obs <- pow(sd(J),2) / mean(J)

  # 2. Compute fit statistic for replicate data
  # 2.1. Discrepancy measure: mean absolute error
  for (i in 1:length(J)){
    J.rep[i] ~ dpois(f[age[i], year[i]])              # Generate replicate data
    D.rep[i] <- J.rep[i] * log(J.rep[i]/J.exp[i]) - (J.rep[i] - J.exp[i])
  }
  Dd.rep <- sum(D.rep)
  # 2.2. Compute test statistic variance-mean ratio
  Tvm.rep <- pow(sd(J.rep),2) / mean(J.rep)

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

  # Assessing the fit of the capture-recapture model
  # 1. Compute fit statistic for observed data
  # 1.1. Discrepancy measure: Freeman-Tukey statistic
  for (t in 1:(n.occasions-1)){
    for (j in 1:n.occasions){
      marr.j.exp[t,j] <- pr.j[t,j] * rel.j[t]         # Expected values
      marr.a.exp[t,j] <- pr.a[t,j] * rel.a[t]         # Expected values
      Dcjs.obs[t,j] <- pow(pow(marr.j[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.obs[t+n.occasions-1,j] <- pow(pow(marr.a[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.obs <- sum(Dcjs.obs)

  # 2. Compute fit statistic for replicate data
  # 2.1. Discrepancy measure: Freeman-Tukey statistic
  for (t in 1:(n.occasions-1)){
    marr.j.rep[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t]) # Generate replicate data
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
inits <- function(){list(mean.p = runif(1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma", "Dmape.obs", "Dmape.rep",
    "Tturn.obs", "Tturn.rep", "Dd.obs", "Dd.rep", "Tvm.obs", "Tvm.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000

# Call JAGS (ART 1 min), check convergence and summarize posteriors
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1)
print(out1, 3)
#               mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# mean.sj      0.301  0.012    0.278    0.300    0.324    FALSE 1 1.000  7500
# mean.sa      0.541  0.013    0.516    0.541    0.566    FALSE 1 1.000  5682
# mean.p       0.606  0.019    0.568    0.606    0.643    FALSE 1 1.000  7307
# mean.f[1]    2.672  0.077    2.524    2.673    2.822    FALSE 1 1.000  7500
# mean.f[2]    3.676  0.087    3.507    3.678    3.846    FALSE 1 1.000  7500
# N[1,1]      44.689 30.084    2.000   40.000  106.000    FALSE 1 1.002   872
# N[2,1]      60.384 27.928    6.000   63.000  105.000    FALSE 1 1.002   773
# N[1,2]      54.207  7.187   40.000   54.000   69.000    FALSE 1 1.001  1946
# ... [output truncated] ...
# N[1,20]     71.730  8.106   56.000   72.000   88.000    FALSE 1 1.001  1384
# N[2,20]     80.911  7.218   67.000   81.000   95.000    FALSE 1 1.000  7500
# sigma       14.149  3.413    8.660   13.713   21.872    FALSE 1 1.000  7500
# Dmape.obs    1.753  0.316    1.196    1.733    2.430    FALSE 1 1.000  7500
# Dmape.rep    1.836  0.610    0.958    1.733    3.291    FALSE 1 1.000  7500
# Tturn.obs   12.000  0.000   12.000   12.000   12.000    FALSE 1   NaN     1
# Tturn.rep   11.384  1.913    8.000   11.000   15.000    FALSE 1 1.000  7500
# Dd.obs     526.625  1.018  525.655  526.301  529.451    FALSE 1 1.000  7500
# Dd.rep     505.500 23.990  458.561  505.476  553.758    FALSE 1 1.000  7500
# Tvm.obs      1.123  0.000    1.123    1.123    1.123    FALSE 1   NaN     1
# Tvm.rep      1.078  0.054    0.976    1.077    1.187    FALSE 1 1.000  7500
# DFT.obs     40.387  1.377   38.306   40.197   43.570    FALSE 1 1.001  2632
# DFT.rep     44.305  4.988   35.164   44.060   54.851    FALSE 1 1.000  7500
# deviance  4291.706  7.552 4276.959 4291.616 4306.746    FALSE 1 1.000  7500

# ~~~~ Plotting function for section 7.2.1 ~~~~

# Do the scatter plot of observed vs simulated discrepancy measures
# jagsout: a jagsUI output object
# obs and rep: character, name of the distrapancy measure to plot
# showP: if TRUE, the Bayesian p-value will be displayed on the plot

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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Do density plots, the denstrip package must be installed ~~~~

# jagsname: character, the name of the jagsUI output object
# param: character, the name of the parameter as in sims.list
# param2: character, the name of the parameter as in summary, if different to param
# ...other arguments as for denstrip.

plotDS <- function(jagsname, param, param2, at, twd=c(3,1.5,1.5),
    tlen=c(2,2,2), width=1/3, colmax="black"){
  if(missing(param2))
    param2 <- param
  ticks <- eval(parse(text=paste0(jagsname, "$summary['", param2, "',c(1,3,7)]")))
  denstrip::denstrip(x=eval(parse(text=paste0(jagsname, "$sims.list$", param))),
      at=at, horiz=FALSE, ticks=ticks, twd=twd, tlen=tlen, width=width, colmax=colmax)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Code for Fig. 7.1 ~~~~
library(scales)
co <- viridis_pal(option='E')(20)[5]

op <- par(mfrow=c(3, 2), las=1)
plotGOF(out1, "Dmape.obs", "Dmape.rep", main="State-space model", col=alpha(co, 0.3))

hist(out1$sims.list$Tturn.rep, nclass=50, col=alpha(co, 0.3),
  xlab="Number of switches", main=NA, axes=FALSE)
abline(v=out1$mean$Tturn.obs, col="red")
axis(1); axis(2)

plotGOF(out1, "Dd.obs", "Dd.rep", main="Poisson regression model", col=alpha(co, 0.3))

hist(out1$sims.list$Tvm.rep, nclass=50, col=alpha(co, 0.3),
    xlab="variance/mean ", main=NA, axes=FALSE)
abline(v=out1$mean$Tvm.obs, col="red")
axis(1); axis(2)

plotGOF(out1, "DFT.obs", "DFT.rep", main="Cormack-Jolly-Seber model", col=alpha(co, 0.3))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 7.2.1.3 Sensitivity of posterior predictive checks to diagnose
#         mis-specified IPMs
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# ~~~~ code for the two mis-specified models ~~~~
# Load data
library(IPMbook)
data(woodchat5)
marr <- marrayAge(woodchat5$ch, woodchat5$age)

# Bundle data
jags.data <- with(woodchat5, list(marr.j=marr[,,1], marr.a=marr[,,2],
    n.occasions=dim(marr)[2], rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
    J=repro[,1], year=repro[,2], age=repro[,3], C=count, pNinit=dUnif(1, 300)))

# Write JAGS model file
cat(file="model2.txt", "
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

  for (t in 1:n.occasions){
    f[1,t] <- mean.f
    f[2,t] <- mean.f
  }

  sigma ~ dunif(0.5, 1000)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
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

  # Productivity data (Poisson regression model)
  for (i in 1:length(J)){
    J[i] ~ dpois(f[age[i], year[i]])
  }

  # Assessing the fit of the Poisson regression model
  # 1. Compute fit statistic for observed data
  # Discrepancy measure: deviance
  for (i in 1:length(J)){
    J.exp[i] <- f[age[i],year[i]]                          # Expected data
    D.obs[i] <- J[i] * log(J[i]/J.exp[i]) - (J[i] - J.exp[i])
  }
  Dd.obs <- sum(D.obs)

  # 2. Compute fit statistic for replicate data
  # Discrepancy measure: mean absolute error
  for (i in 1:length(J)){
    J.rep[i] ~ dpois(f[age[i], year[i]])                   # Generate replicate data
    D.rep[i] <- J.rep[i] * log(J.rep[i]/J.exp[i]) - (J.rep[i] - J.exp[i])
  }
  Dd.rep <- sum(D.rep)

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
inits <- function(){list(mean.p = runif(1))}

# Parameters monitored
parameters <- c("mean.s", "mean.p", "mean.f", "N", "sigma",
    "Dmape.obs", "Dmape.rep", "Dd.obs", "Dd.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000

# Call JAGS (ART 1 min), check convergence and summarize posteriors
out2 <- jags(jags.data, inits, parameters, "model2.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out2)
print(out2, 3)

# IPM3: no reproduction of 1y shrikes
# Bundle data
jags.data <- with(woodchat5, list(marr.j=marr[,,1], marr.a=marr[,,2],
    n.occasions=dim(marr)[2], rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
    J=repro[repro[,3]==2,1], year=repro[repro[,3]==2,2], C=count, pNinit=dUnif(1, 300)))

# Write JAGS model file
cat(file="model3.txt", "
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

  for (t in 1:n.occasions){
    f[t] <- mean.f
  }

  sigma ~ dunif(0.5, 1000)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois(N[2,t] * f[t] / 2 * sj[t])
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

  # Productivity data (Poisson regression model)
  for (i in 1:length(J)){
    J[i] ~ dpois(f[year[i]])
  }

  # Assessing the fit of the Poisson regression model
  # 1. Compute fit statistic for observed data
  # Discrepancy measure: deviance
  for (i in 1:length(J)){
    J.exp[i] <- f[year[i]]                                 # Expected data
    D.obs[i] <- J[i] * log(J[i]/J.exp[i]) - (J[i] - J.exp[i])
  }
  Dd.obs <- sum(D.obs)

  # 2. Compute fit statistic for replicate data
  # Discrepancy measure: mean absolute error
  for (i in 1:length(J)){
    J.rep[i] ~ dpois(f[year[i]])                           # Generate replicate data
    D.rep[i] <- J.rep[i] * log(J.rep[i]/J.exp[i]) - (J.rep[i] - J.exp[i])
  }
  Dd.rep <- sum(D.rep)

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
inits <- function(){list(mean.p = runif(1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma",
    "Dmape.obs", "Dmape.rep", "Dd.obs", "Dd.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000

# Call JAGS (ART 1 min), check convergence and summarize posteriors
out3 <- jags(jags.data, inits, parameters, "model3.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out3)
print(out3, 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Analysis of the five times larger data set ~~~~
# Load data
library(IPMbook)
data(woodchat7)
marr <- marrayAge(woodchat7$ch, woodchat7$age)

# Bundle data
jags.data <- with(woodchat7, list(marr.j=marr[,,1], marr.a=marr[,,2],
    n.occasions=dim(marr)[2], rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
    J=repro[,1], year=repro[,2], age=repro[,3], C=count, pNinit=dUnif(100, 1000)))

# Fit IPM1
# Initial values
inits <- function(){list(mean.p = runif(1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma",
    "Dmape.obs", "Dmape.rep", "Tturn.obs", "Tturn.rep", "Dd.obs", "Dd.rep",
    "Tvm.obs", "Tvm.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000

# Call JAGS (ART 8 min), check convergence and summarize posteriors
out4 <- jags(jags.data, inits, parameters, "model1.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out4)
print(out4, 3)

# Fit IPM2
# Initial values
inits <- function(){list(mean.p = runif(1))}

# Parameters monitored
parameters <- c("mean.s", "mean.p", "mean.f", "N", "sigma",
    "Dmape.obs", "Dmape.rep", "Dd.obs", "Dd.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000

# Call JAGS (ART 8 min), check convergence and summarize posteriors
out5 <- jags(jags.data, inits, parameters, "model2.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out5)
print(out5, 3)

# Fit IPM3
# Bundle data
marr <- marrayAge(woodchat7$ch, woodchat7$age)
jags.data <- with(woodchat7, list(marr.j=marr[,,1], marr.a=marr[,,2],
    n.occasions=dim(marr)[2], rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
    J=repro[repro[,3]==2,1], year=repro[repro[,3]==2,2], C=count,
    pNinit=dUnif(100, 1000)))

# Initial values
inits <- function(){list(mean.p = runif(1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma",
    "Dmape.obs", "Dmape.rep", "Dd.obs", "Dd.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000

# Call JAGS (ART 8 min), check convergence and summarize posteriors
out6 <- jags(jags.data, inits, parameters, "model3.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out6)
print(out6, 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Separate analyses of capture-recapture and productivity data ~~~~
# CJS and Poisson regression models
# Write JAGS model file
cat(file = "model4.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  for (t in 1:n.occasions){
    f[1,t] <- mean.f[1]
    f[2,t] <- mean.f[2]
  }

  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f[1] ~ dunif(0, 10)
  mean.f[2] ~ dunif(0, 10)

  # Productivity data (Poisson regression model)
  for (i in 1:length(J)){
    J[i] ~ dpois(f[age[i], year[i]])
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
}
")

# Small data set
# Bundle data
marr <- marrayAge(woodchat5$ch, woodchat5$age)
jags.data <- with(woodchat5, list(marr.j=marr[,,1], marr.a=marr[,,2],
    n.occasions=dim(marr)[2], rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
    J=repro[,1], year=repro[,2], age=repro[,3]))

# Initial values
inits <- function(){list(mean.sj = runif(1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f")

# MCMC settings
ni <- 12000; nb <- 2000; nt <- 4; nc <- 3; na <- 2000

# Call JAGS from R (ART <1 min) and check convergence
out7 <- jags(jags.data, inits, parameters, "model4.txt",
    n.iter = ni, n.burnin = nb, n.chains = nc, n.thin = nt, n.adapt = na, parallel = TRUE)

# Large data set
data(woodchat7)
marr <- marrayAge(woodchat7$ch, woodchat7$age)
jags.data <- with(woodchat7, list(marr.j=marr[,,1], marr.a=marr[,,2],
    n.occasions=dim(marr)[2], rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
    J=repro[,1], year=repro[,2], age=repro[,3]))

# Initial values
inits <- function(){list(mean.sj=runif(1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f")

# MCMC settings
ni <- 12000; nb <- 2000; nt <- 4; nc <- 3; na <- 2000

# Call JAGS from R (ART 3 min) and check convergence
out8 <- jags(jags.data, inits, parameters, "model4.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

save(out1, out2, out3, out4, out5, out6, out7, out8, file="ResultsChapter7.2.1.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for Figure 7.2 ~~~~
load("ResultsChapter7.2.1.Rdata")
library(scales)
library(plotrix)
co <- viridis_pal(option='E')(20)[c(5, 11, 16)]

op <- par(las=1, mar=c(3,7,3,1), "mfrow")
layout(matrix(1:9, 3, 3, byrow=TRUE), widths=c(1.1, 1, 1), heights=c(1, 1, 1.1), TRUE)

plotGOF(out1, "Dmape.obs", "Dmape.rep", main="State-space model", col=alpha(co[1], 0.3))
mtext(expression(bold(IPM[1])), side=2, las=0, line=5, font=2)
corner.label('A', font=2, cex=1.25)

par(mar=c(3,4,3,1))
plotGOF(out1, "Dd.obs", "Dd.rep", main="Poisson regression model", col=alpha(co[1], 0.3))
corner.label('B', font=2, cex=1.25)

plotGOF(out1, "DFT.obs", "DFT.rep", main="Cormack-Jolly-Seber model", col=alpha(co[1], 0.3))
corner.label('C', font=2, cex=1.25)

par(mar=c(3,7,3,1))
plotGOF(out2, "Dmape.obs", "Dmape.rep", main="State-space model", col=alpha(co[2], 0.3))
mtext(expression(bold(IPM[2])), side=2, las=0, line=5, font=2)
corner.label('D', font=2, cex=1.25)

par(mar=c(3,4,3,1))
plotGOF(out2, "Dd.obs", "Dd.rep", main="Poisson regression model", col=alpha(co[2], 0.3))
corner.label('E', font=2, cex=1.25)

plotGOF(out2, "DFT.obs", "DFT.rep", main="Cormack-Jolly-Seber model", col=alpha(co[2], 0.3))
corner.label('F', font=2, cex=1.25)

par(mar=c(5,7,3,1))
plotGOF(out3, "Dmape.obs", "Dmape.rep", main="State-space model", col=alpha(co[3], 0.3))
mtext(expression(bold(IPM[2])), side=2, las=0, line=5, font=2)
corner.label('G', font=2, cex=1.25)

par(mar=c(5,4,3,1))
plotGOF(out3, "Dd.obs", "Dd.rep", main="Poisson regression model", col=alpha(co[3], 0.3))
corner.label('H', font=2, cex=1.25)

plotGOF(out3, "DFT.obs", "DFT.rep", main="Cormack-Jolly-Seber model", col=alpha(co[3], 0.3))
corner.label('I', font=2, cex=1.25)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for Fig. 7.3 ~~~~

load("ResultsChapter7.2.1.Rdata")

library(scales)
library(plotrix)
co <- viridis_pal(option='E')(20)[c(5, 11, 16)]

op <- par(las=1, mar=c(3,7,3,1), "mfrow")
layout(matrix(1:9, 3, 3, byrow=TRUE), widths=c(1.1, 1, 1), heights=c(1, 1, 1.1), TRUE)
plotGOF(out4, "Dmape.obs", "Dmape.rep", main="State-space model", col=alpha(co[1], 0.3))
mtext(expression(bold(IPM[1])), side=2, las=0, line=5, font=2)
corner.label('A', font=2, cex=1.25)

par(mar=c(3,4,3,1))
plotGOF(out4, "Dd.obs", "Dd.rep", main="Poisson regression model", col=alpha(co[1], 0.3))
corner.label('B', font=2, cex=1.25)

plotGOF(out4, "DFT.obs", "DFT.rep", main="Cormack-Jolly-Seber model", col=alpha(co[1], 0.3))
corner.label('C', font=2, cex=1.25)

par(mar=c(3,7,3,1))
plotGOF(out5, "Dmape.obs", "Dmape.rep", main="State-space model", col=alpha(co[2], 0.3))
mtext(expression(bold(IPM[2])), side=2, las=0, line=5, font=2)
corner.label('D', font=2, cex=1.25)

par(mar=c(3,4,3,1))
plotGOF(out5, "Dd.obs", "Dd.rep", main="Poisson regression model", col=alpha(co[2], 0.3))
corner.label('E', font=2, cex=1.25)

plotGOF(out5, "DFT.obs", "DFT.rep", main="Cormack-Jolly-Seber model", col=alpha(co[2], 0.3))
corner.label('F', font=2, cex=1.25)

par(mar=c(5,7,3,1))
plotGOF(out6, "Dmape.obs", "Dmape.rep", main="State-space model", col=alpha(co[3], 0.3))
mtext(expression(bold(IPM[2])), side=2, las=0, line=5, font=2)
corner.label('G', font=2, cex=1.25)

par(mar=c(5,4,3,1))
plotGOF(out6, "Dd.obs", "Dd.rep", main="Poisson regression model", col=alpha(co[3], 0.3))
corner.label('H', font=2, cex=1.25)

plotGOF(out6, "DFT.obs", "DFT.rep", main="Cormack-Jolly-Seber model", col=alpha(co[3], 0.3))
corner.label('I', font=2, cex=1.25)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for Figure 7.4 ~~~~

load("ResultsChapter7.2.1.Rdata")

library(denstrip)
library(scales)
co <- viridis_pal(option='E')(20)[c(5, 11, 16)]
labs <- c(expression(IPM[1]), expression(IPM[2]), expression(IPM[3]),
    expression(CJS[ ]), expression(IPM[1]), expression(IPM[2]),
    expression(IPM[3]), expression(CJS[ ]))

op <- par(las=1, mar=c(5, 5, 2, 2), "mfrow")
layout(matrix(1:4, 2, 2, byrow=TRUE), widths=c(1.25, 1.25), heights=c(1, 1), TRUE)

# Top row
plot(0, ylim=c(0.25, 0.45), xlim=c(0.7,9.3), axes=F, pch=NA, xlab=NA,
    ylab="Juvenile survival")
plotDS("out1", "mean.sj", at=1, colmax=co[1])
plotDS("out2", "mean.s", at=2, colmax=co[2])
plotDS("out3", "mean.sj", at=3, colmax=co[3])
plotDS("out7", "mean.sj", at=4)
plotDS("out4", "mean.sj", at=6, colmax=co[1])
plotDS("out5", "mean.s", at=7, colmax=co[2])
plotDS("out6", "mean.sj", at=8, colmax=co[3])
plotDS("out8", "mean.sj", at=9)
axis(2, las=1)
axis(1, at=c(1:4, 6:9), labels=labs)
mtext("Small sample size", side=1, at=2.5, line=3)
mtext("Large sample size", side=1, at=7.5, line=3)
corner.label('A', font=2, xoff=0.25)

plot(0, ylim=c(0.4, 0.6), xlim=c(0.7,9.3), axes=F, pch=NA, xlab=NA,
    ylab="Adult survival")
plotDS("out1", "mean.sa", at=1, colmax=co[1])
plotDS("out2", "mean.s", at=2, colmax=co[2])
plotDS("out3", "mean.sa", at=3, colmax=co[3])
plotDS("out7", "mean.sa", at=4)
plotDS("out4", "mean.sa", at=6, colmax=co[1])
plotDS("out5", "mean.s", at=7, colmax=co[2])
plotDS("out6", "mean.sa", at=8, colmax=co[3])
plotDS("out8", "mean.sa", at=9)
axis(2, las=1)
axis(1, at=c(1:4, 6:9), labels=labs)
mtext("Small sample size", side=1, at=2.5, line=3)
mtext("Large sample size", side=1, at=7.5, line=3)
corner.label('B', font=2, xoff=0.25)

# Second row
labs[c(4,8)] <- expression(Pois[ ])
plot(0, ylim=c(2.4, 3.3), xlim=c(0.7,9.3), axes=F, pch=NA, xlab=NA,
    ylab="First year productivity")
plotDS("out1", "mean.f[,1]", "mean.f[1]", at=1, colmax=co[1])
plotDS("out2", "mean.f", at=2, colmax=co[2])
plotDS("out7", "mean.f[,1]", "mean.f[1]", at=4)
plotDS("out4", "mean.f[,1]", "mean.f[1]", at=6, colmax=co[1])
plotDS("out5", "mean.f", at=7, colmax=co[2])
plotDS("out8", "mean.f[,1]", "mean.f[1]", at=9)
axis(2, las=1)
axis(1, at=c(1:4, 6:9), labels=labs)
mtext("Small sample size", side=1, at=2.5, line=3)
mtext("Large sample size", side=1, at=7.5, line=3)
corner.label('C', font=2, xoff=0.25)

plot(0, ylim=c(2.95, 4.05), xlim=c(0.7,9.3), axes=F, pch=NA, xlab=NA,
    ylab="Adult productivity")
plotDS("out1", "mean.f[,2]", "mean.f[2]", at=1, colmax=co[1])
plotDS("out2", "mean.f", at=2, colmax=co[2])
plotDS("out3", "mean.f", at=3, colmax=co[3])
plotDS("out7", "mean.f[,2]", "mean.f[2]", at=4)
plotDS("out4", "mean.f[,2]", "mean.f[2]", at=6, colmax=co[1])
plotDS("out5", "mean.f", at=7, colmax=co[2])
plotDS("out6", "mean.f", at=8, colmax=co[3])
plotDS("out8", "mean.f[,2]", "mean.f[2]", at=9)
axis(2, las=1)
axis(1, at=c(1:4, 6:9), labels=labs)
mtext("Small sample size", side=1, at=2.5, line=3)
mtext("Large sample size", side=1, at=7.5, line=3)
corner.label('D', font=2, xoff=0.25)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for Figure 7.5 ~~~~

load("ResultsChapter7.2.1.Rdata")

library(IPMbook)
data(woodchat5)
data(woodchat7)

library(scales)
co <- viridis_pal(option='E')(20)[c(5, 11, 16)]

# Produce plot with counts against estimated population size from the 3 models
Ntot <- NtotL <- array(NA, dim=c(dim(out1$sims.list$N)[1], 3, dim(out1$sims.list$N)[3]))
for (t in 1:length(woodchat5$count)){
  Ntot[,1,t] <- out1$sims.list$N[,1,t] + out1$sims.list$N[,2,t]
  Ntot[,2,t] <- out2$sims.list$N[,1,t] + out2$sims.list$N[,2,t]
  Ntot[,3,t] <- out3$sims.list$N[,1,t] + out3$sims.list$N[,2,t]

  NtotL[,1,t] <- out4$sims.list$N[,1,t] + out4$sims.list$N[,2,t]
  NtotL[,2,t] <- out5$sims.list$N[,1,t] + out5$sims.list$N[,2,t]
  NtotL[,3,t] <- out6$sims.list$N[,1,t] + out6$sims.list$N[,2,t]
}
qu <- function(x) quantile(x, c(0.025, 0.975))
T <- length(woodchat5$count)
d <- 0.2

op <- par(las=1, mar=c(3,5,3,1), "mfrow")
layout(matrix(1:2, 2, 1, byrow=TRUE), widths=1.6, heights=c(1, 1), TRUE)

plot(y=apply(Ntot[,1,], 2, mean), x=(1:T)-d, type="b", ylim=c(80, 200),
    ylab="Number", xlab=NA, las=1, pch=16, axes=FALSE, col=co[1], main="Small sample size")
axis(2)
axis(1, at=1:20, tcl=-0.25, labels=NA)
axis(1, at=c(5, 10, 15, 20), tcl=-0.5, labels=NA)

segments((1:T)-d, apply(Ntot[,1,], 2, qu)[1,], (1:T)-d, apply(Ntot[,1,], 2, qu)[2,], col=co[1])
points(y=apply(Ntot[,2,], 2, mean), x=1:T, type="b", pch=16, col=co[2])
segments(1:T, apply(Ntot[,2,], 2, qu)[1,], 1:T, apply(Ntot[,2,], 2, qu)[2,], col=co[2])
points(y=apply(Ntot[,3,], 2, mean), x=(1:T)+d, type="b", pch=16, col=co[3])
segments((1:T)+d, apply(Ntot[,3,], 2, qu)[1,], (1:T)+d, apply(Ntot[,3,], 2, qu)[2,], col=co[3])
points(woodchat5$count, type="b", pch=1, lty=3)
legend("topleft", legend=c("Observed counts", expression("Estimates from IPM"[1]),
    expression("Estimates from IPM"[2]), expression("Estimates from IPM"[3])),
    pch=c(1, rep(16, 3)), col=c("black", co), lty=c(3, 1, 1, 1), bty="n")

par(mar=c(5,5,1,1))
plot(y=apply(NtotL[,1,], 2, mean), x=(1:T)-d, type="b", ylim=c(0, 1600), ylab="Number",
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
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
