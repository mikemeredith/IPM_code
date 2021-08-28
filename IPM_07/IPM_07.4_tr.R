# Schaub & Kéry (2022) Integrated Population Models
# Chapter 7 : Assessment of integrated population models
# ------------------------------------------------------
# Code from final MS.

# Run time for test script 2 mins, full run 3 hrs

# 7.4 Effects of a mis-specified observation model
# ================================================

# ~~~ Code for the simulations of section 7.4 ~~~

library(jagsUI)
library(IPMbook)

# Analysing models

# IPM1: with all 3 data sets
cat(file="model8bis.txt", "
model {
  # Priors and constraints
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

  sigma.obs ~ dunif(0.5, upper.sigma.obs)
  tau.obs <- pow(sigma.obs, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois((N[1,t] + N[2,t]) * f[t] / 2 * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)
  }

  # Assessing the fit of the state-space model
  for (t in 1:n.occasions){
    C.exp[t] <- N[1,t] + N[2,t]                    # Expected counts
    Dssm.obs[t] <- abs((C[t] - C.exp[t]) / C[t])   # Discrepancy measure

    C.rep[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)     # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t]) # Discrepancy measure
  }
  Dmape.obs <- sum(Dssm.obs)
  Dmape.rep <- sum(Dssm.rep)


  # Productivity data (Poisson regression model)
  for (t in 1:n.occasions){
    J[t] ~ dpois(B[t] * f[t])

    J.exp[t] <- B[t] * f[t]              # Expected data
    D.obs[t] <- J[t] * log(J[t]/J.exp[t]) - (J[t] - J.exp[t])

    J.rep[t] ~ dpois(B[t] * f[t])
    D.rep[t] <- J.rep[t] * log(J.rep[t]/J.exp[t]) - (J.rep[t] - J.exp[t])
  }
  Dd.obs <- sum(D.obs)
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
    q[t] <- 1-p[t]   # Probability of non-recapture
    pr.j[t,t] <- sj[t]*p[t]
    pr.a[t,t] <- sa[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t]*prod(sa[(t+1):j])*prod(q[t:(j-1)])*p[j]
      pr.a[t,j] <- prod(sa[t:j])*prod(q[t:(j-1)])*p[j]
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

  # Assessing the fit of the capture-recapture model (Freeman-Tukey)
  for (t in 1:(n.occasions-1)){
    marr.j.rep[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])   # Generate replicate data
    marr.a.rep[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    for (j in 1:n.occasions){
      marr.j.exp[t,j] <- pr.j[t,j] * rel.j[t]  # Expected values
      marr.a.exp[t,j] <- pr.a[t,j] * rel.a[t]  # Expected values
      Dcjs.obs[t,j] <- pow(pow(marr.j[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.obs[t+n.occasions-1,j] <- pow(pow(marr.a[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)

      Dcjs.rep[t,j] <- pow(pow(marr.j.rep[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.rep[t+n.occasions-1,j] <- pow(pow(marr.a.rep[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.obs <- sum(Dcjs.obs)
  DFT.rep <- sum(Dcjs.rep)

  # Derived parameters
  # Mean population growth rate (geometric mean)
  geom.rate <- pow(Ntot[n.occasions] / Ntot[1], 1 / (n.occasions-1))
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")



# IPM2: with counts and CMR only
cat(file="model9bis.txt", "
model {
  # Priors and constraints
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

  sigma.obs ~ dunif(0.5, upper.sigma.obs)
  tau.obs <- pow(sigma.obs, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois((N[1,t] + N[2,t]) * f[t] / 2 * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)
  }

  # Assessing the fit of the state-space model
  for (t in 1:n.occasions){
    C.exp[t] <- N[1,t] + N[2,t]                    # Expected counts
    Dssm.obs[t] <- abs((C[t] - C.exp[t]) / C[t])   # Discrepancy measure

    C.rep[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)     # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t]) # Discrepancy measure
  }
  Dmape.obs <- sum(Dssm.obs)
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
    q[t] <- 1-p[t]   # Probability of non-recapture
    pr.j[t,t] <- sj[t]*p[t]
    pr.a[t,t] <- sa[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t]*prod(sa[(t+1):j])*prod(q[t:(j-1)])*p[j]
      pr.a[t,j] <- prod(sa[t:j])*prod(q[t:(j-1)])*p[j]
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

  # Assessing the fit of the capture-recapture model (Freeman-Tukey)
  for (t in 1:(n.occasions-1)){
    marr.j.rep[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])   # Generate replicate data
    marr.a.rep[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    for (j in 1:n.occasions){
      marr.j.exp[t,j] <- pr.j[t,j] * rel.j[t]  # Expected values
      marr.a.exp[t,j] <- pr.a[t,j] * rel.a[t]  # Expected values
      Dcjs.obs[t,j] <- pow(pow(marr.j[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.obs[t+n.occasions-1,j] <- pow(pow(marr.a[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)

      Dcjs.rep[t,j] <- pow(pow(marr.j.rep[t,j], 0.5) - pow(marr.j.exp[t,j], 0.5), 2)
      Dcjs.rep[t+n.occasions-1,j] <- pow(pow(marr.a.rep[t,j], 0.5) - pow(marr.a.exp[t,j], 0.5), 2)
    } #j
  } #t
  DFT.obs <- sum(Dcjs.obs)
  DFT.rep <- sum(Dcjs.rep)

  # Derived parameters
  # Mean population growth rate (geometric mean)
  geom.rate <- pow(Ntot[n.occasions] / Ntot[1], 1 / (n.occasions-1))
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")


# Simulations

# Number of simulations
# nsim <- 600
nsim <- 6  # ~~~ for testing

# Number of years
T <- 10

# Observation error for the population survey
sigma <- 10
p.const <- 0.6
p.trend <- seq(0.3, 0.7, length.out=T)

# Capture and recapture probabilities
cap <- 0.4                        # initial capture probability
recap <- 0.6                      # recapture probability

# Probability to find a brood whose reproductive ouput is recorded
pprod <- 0.5

# Initial population size per age class
Ni <- c(50, 50)

# Demographic rates
# Survival
phi <- c(0.3, 0.55)

# Productivity
f <- 3.1

# Define matrices to store results
res3 <- array(NA, dim=c(43, 11, nsim, 4))
res2 <- array(NA, dim=c(41, 11, nsim, 4))
bP3 <- array(NA, dim=c(nsim, 3, 4))
bP2 <- array(NA, dim=c(nsim, 2, 4))
Nsize <- matrix(NA, ncol=T, nrow=nsim)

# Initial values
inits <- function(){
  N <- matrix(NA, nrow=2, ncol=T)
  N[1,] <- round(runif(T, 40, 60))
  N[2,] <- round(runif(T, 40, 60))
  ini <- list(mean.sj=runif(1, 0, 0.5), N=N)
  return(ini)
}

# Parameters monitored
parameters1 <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "geom.rate",
    "Ntot", "sigma.obs", "Dmape.obs", "Dmape.rep", "DFT.obs", "DFT.rep",
    "Dd.obs", "Dd.rep")
parameters2 <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "geom.rate",
    "Ntot", "sigma.obs", "Dmape.obs", "Dmape.rep", "DFT.obs", "DFT.rep")

# MCMC settings
ni <- 4000; nb <- 2000; nc <- 3; nt <- 4; na <- 2000


# Start simulations
system.time(
for (s in 1:nsim){

  set.seed(s)

  # Create populations
  pop1 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
  pop2 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
  pop3 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)

  # Create capture histories & m-arrays from population 1
  ch <- simCapHist(state=pop1$state, cap=cap, recap=recap, maxAge=2)
  marr <- marrayAge(ch$ch, ch$age)

  # Create productivity data from population 2
  pro <- simProd(reprod=pop2$reprod, pInclude=pprod)
  # Aggregate productivity data to make the model run faster
  J <- pro$prod.agg[,1]
  B <- pro$prod.agg[,2]

  # Create the population survey data from population 3
  # a: normal distibution
  count1 <- simCountNorm(N=pop3$totB, sigma=sigma)$count
  # b: binomial distribution with constant detection
  count2 <- simCountBin(N=pop3$totB, pDetect=p.const)$count
  # c: binomial distribution with temporally variable detection
  p.t <- runif(T, 0.4, 0.8)
  count3 <- simCountBin(N=pop3$totB, pDetect=p.t)$count
  # d: binomial distribution with detection that has a trend
  count4 <- simCountBin(N=pop3$totB, pDetect=p.trend)$count

  Nsize[s,] <- pop3$totB   # Store true population size

  # Bundle data
  jagsData <- vector('list', 4)
  # a: use ‘Normal’ count
  jagsData[[1]] <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
      C=count1, pNinit=dUnif(1, 100), upper.sigma.obs=100, J=J, B=B)

  # b: use ‘binomial’ count with constant detection
  jagsData[[2]] <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
      C=count2, pNinit=dUnif(1, 100), upper.sigma.obs=100, J=J, B=B)

  # c: use ‘binomial’ count with variable detection
  jagsData[[3]] <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
      C=count3, pNinit=dUnif(1, 100), upper.sigma.obs=100, J=J, B=B)

  # d: use ‘binomial’ count with detection that has a trend
  jagsData[[4]] <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]),
      C=count4, pNinit=dUnif(1, 100), upper.sigma.obs=100, J=J, B=B)

  # Call JAGS from R (jagsUI)
  for (i in 1:4){
    # 3 data sets
    out <- try(jags(jagsData[[i]], inits, parameters1, "model8bis.txt",
        n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
    if(!inherits(out, "try-error")){
      res3[,,s, i] <- out$summary
      bP3[s,1,i] <- mean(out$sims.list$Dmape.rep > out$sims.list$Dmape.obs)
      bP3[s,2,i] <- mean(out$sims.list$DFT.rep > out$sims.list$DFT.obs)
      bP3[s,3,i] <- mean(out$sims.list$Dd.rep > out$sims.list$Dd.obs)
    } # if

    # 2 data sets (without J and B)
    out <- try(jags(jagsData[[i]][1:8], inits, parameters2, "model9bis.txt",
        n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
    if(!inherits(out, "try-error")){
      res2[,,s,i] <- out$summary
      bP2[s,1,i] <- mean(out$sims.list$Dmape.rep > out$sims.list$Dmape.obs)
      bP2[s,2,i] <- mean(out$sims.list$DFT.rep > out$sims.list$DFT.obs)
    } # if
  } # i
  print(s)
} ) # s


save(res3, res2, phi, f, sigma, cap, recap, p.const, p.trend, Nsize, bP3, bP2,
    file="ResultsChapter7.4.Rdata")

# Produce plots
load('ResultsChapter7.4.Rdata')
str(res3)  # nodes x stats x sims x data (data = 4 scenarios for count data)
str(res2)  # same but fewer nodes

# Select only converged simulations, based on Rhat for params:
tocheck <- c(1, 2, 4, 25:35)
# Take a look at last simulation to see what these are:
out$summary[tocheck, 8]

r.crit <- 1.05
# the MS uses 1.1 for Fig 7.16, 1.05 for the others.
rhats3 <- res3[tocheck,8,,] < r.crit
good3 <- apply(rhats3, 2, all)
rhats2 <- res2[tocheck,8,,] < r.crit
good2 <- apply(rhats2, 2, all)
incl <- which(good3 & good2)
if(length(incl) > 500)  # ~~~ will be fewer when testing
  incl <- incl[1:500]


# Boxplots for 5 parameters plus one blank plot
parstoplot <- c(1, 2, NA, 4, 25, 36)
# Take a look at last simulation to see what these are:
out$summary[parstoplot,1:2]

# True values are all scalar
( true <- c(phi, NA, f, phi[1] * f / 2 + phi[2], 10) )

toplot3 <- res3[parstoplot, , incl,]
toplot2 <- res2[parstoplot, , incl,]

# Get relative bias
# est is parameters x stats x sims x data, true is vector of length parameters
rbias <- function(est, true){
  means <- est[,1,,]  # 6 x 50 x 4
  errs <- sweep(means, 1, true, "-")
  return(sweep(errs, 1, true, "/"))
}

all_1 <- abind::abind(rbias(toplot3, true), rbias(toplot2, true), along=3)

# Sort out colours, names, etc
library(scales)
co <- viridis_pal(option='E')(20)[c(2, 7, 12, 17)]

xnames <- c("N", expression('p'[.]), expression('p'[t]), expression('p'[T]),
    "N", expression('p'[.]), expression('p'[t]), expression('p'[T]))
ylabs <- c(expression('Relative bias ('*italic(s)[j]*')'),
           expression('Relative bias ('*italic(s)[a]*')'),
           NA,
           expression('Relative bias ('*italic(f)*')'),
           expression('Relative bias ('*lambda*')'),
           expression('Relative bias ('*sigma*')'))


# Fig. 7.16
op <- par(mfrow=c(3, 2), las=1, mar=c(4, 5, 0, 1))
for (i in 1:6){
  if(i == 3){
    # plot.new()
    next
  }
  boxplot(all_1[i,,], ylab=ylabs[i], outline=FALSE, xlab=NA,
      col=c(co, rep("white", 4)), border=c(rep("black", 4), co),
      lwd=c(rep(1,4), rep(2,4)), axes=FALSE, boxwex=0.6)
  if(i > 4){
    axis(1, at=1:8, labels=xnames)
  }
  else{
    axis(1, at=1:8, labels=NA)
  }
  axis(2)
} # i
par(op)

############

# Fig. 7.17

load('ResultsChapter7.4.Rdata')

rmse <- function(est, true){
  means <- est[,1,,]  # 6 x 50 x 4
  SDs <- est[,2,,]  # 6 x 50 x 4
  errs <- sweep(means, 1, true, "-")
  return(sqrt(errs^2 + SDs^2))
}

all_2 <- abind::abind(rmse(toplot3, true), rmse(toplot2, true), along=3)

ylabs <- c(expression('RMSE ('*italic(s)[j]*')'),
           expression('RMSE ('*italic(s)[a]*')'),
           NA,
           expression('RMSE ('*italic(f)*')'),
           expression('RMSE ('*lambda*')'),
           expression('RMSE ('*sigma*')'))

op <- par(mfrow=c(3, 2), las=1, mar=c(4, 5, 0, 1))
for (i in 1:6){
  if(i == 3)
    next
  boxplot(all_2[i,,], ylab=ylabs[i], outline=FALSE, xlab=NA, col=c(co, rep("white", 4)),
      border=c(rep("black", 4), co), lwd=c(rep(1,4), rep(2,4)), axes=FALSE, boxwex=0.6)
  if(i > 4){
    axis(1, at=1:8, labels=xnames)
  }
  else{
    axis(1, at=1:8, labels=NA)
  }
  axis(2)
} # i
par(op)

##############


# Fig. 7.18

load('ResultsChapter7.4.Rdata')

parstoplot <- 26:35
# Take a look at last simulation to see what these are:
out$summary[parstoplot, 1:2]
# We only need the means
toplot3 <- res3[parstoplot, 1, incl,]
toplot2 <- res2[parstoplot, 1, incl,]

true.N <- apply(Nsize[incl,], 2, mean)
qu <- function(x) quantile(x, p=c(0.025, 0.975))
lower.N <- apply(Nsize[incl,], 2, qu)[1,]
upper.N <- apply(Nsize[incl,], 2, qu)[2,]

library(scales)
co <- viridis_pal(option='E')(20)[c(2, 7, 12, 17)]

ylim <- c(20, 190)
T <- 10
p <- cbind(NA, p.const, 0.6, p.trend)

leftLabels <- c(
    expression(bold('True observation model: N')),
    expression(bold('True observation model: p'[.])),
    expression(bold('True observation model: p'[t])),
    expression(bold('True observation model: p'[T])) )

op <- par(las=1, mfrow=c(4, 2), oma=c(2.5, 2, 2, 0), "mar")
for (i in 1:4){
  par(mar=c(1, 5, 1, 0))
  boxplot(t(toplot3[,,i]), outline=FALSE, col=co[i], axes=FALSE, ylim=ylim,
      boxwex=0.5, xpd=TRUE, ylab="Population size")
  axis(1, at=1:10, tcl=-0.25, labels=NA)
  axis(1, at=c(1, 3, 5, 7, 9), tcl=-0.5, labels=(i==4))
  axis(2)
  title(ylab="Population size", xpd=TRUE)
  points(true.N, pch=16, col="orange")
  segments(1:T, lower.N, 1:T, upper.N, col="orange")
  lines(p[,i] * true.N, lty=2, col="magenta")
  if(i == 1){
    legend('bottomleft', pch=16, col="orange", legend='True population size', bty='n')
    mtext(expression(bold(IPM[1])), side=3, line=1)
  } # if
  if(i == 2){
    legend('bottomleft', lty=2, col="magenta", legend='Expected population size', bty='n')
  } # if
  mtext(leftLabels[i], side=2, line=5, las=3)
  if(i == 4){
    mtext('Time', side=1, line=2.2)
  } # if
  par(mar=c(1, 3, 1, 2))
  boxplot(t(toplot2[,,i]), ylab=NA, outline=FALSE, col="white", border=co[1], lwd=2,
      axes=FALSE, ylim=ylim, xlab=NA, boxwex=0.5)
  axis(1, at=1:10, tcl=-0.25, labels=NA)
  axis(1, at=c(1, 3, 5, 7, 9), tcl=-0.5, labels=(i==4))
  axis(2, labels=NA)
  points(true.N, pch=16, col="orange")
  segments(1:T, lower.N, 1:T, upper.N, col="orange")
  lines(p[,i] * true.N, lty=2, col="magenta")
  if(i == 1){
    mtext(expression(bold(IPM[2])), side=3, line=1)
  } # if
  if(i == 4){
    mtext('Time', side=1, line=2.2)
  } # if
} # i
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
