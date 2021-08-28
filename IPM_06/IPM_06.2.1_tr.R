# Schaub & Kéry (2022) Integrated Population Models
# Chapter 6 : Benefits of integrated population modeling
# ------------------------------------------------------
# Code from final MS.

# Run time testing 40 secs, full run 2.2 hrs


# 6.2 Parameter estimates with increased precision
# ================================================

# 6.2.1 Experiencing the gain in precision in a simple simulation
# ---------------------------------------------------------------

# ~~~~ Code for the simulations ~~~~

library(IPMbook); library(jagsUI)

# 1. Write JAGS code for the analyzing models
# 1.1. IPM
cat(file="model6.txt", "
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

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- mean.f/2 * mean.sj * (N[1,t] + N[2,t])
    N[2,t+1] <- mean.sa * (N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Productivity data (Poisson regression model)
  sJ ~ dpois(nJ * mean.f)

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


# 1.2. CJS model
cat(file="model7.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)

  for (t in 1:(n.occCJS-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  # Define the multinomial likelihood
  for (t in 1:(n.occCJS-1)){
    marr.j[t,1:n.occCJS] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occCJS] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occCJS-1)){
    # Main diagonal
    q[t] <- 1 - p[t]   # Probability of non-recapture
    pr.j[t,t] <- sj[t] * p[t]
    pr.a[t,t] <- sa[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occCJS-1)){
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
  for (t in 1:(n.occCJS-1)){
    pr.j[t,n.occCJS] <- 1-sum(pr.j[t,1:(n.occCJS-1)])
    pr.a[t,n.occCJS] <- 1-sum(pr.a[t,1:(n.occCJS-1)])
  }
}
")


# 1.3. Poisson regression model
cat(file="model8.txt", "
model {
  # Priors and linear models
  mean.f ~ dunif(0, 10)

  # Likelihood
  sJ ~ dpois(nJ * mean.f)
}
")


# 1.4. State-space model
cat(file="model9.txt", "
model {
  # Priors and linear models
  lambda ~ dunif(0, 5)

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Model for the initial population size: uniform priors
  N[1] ~ dunif(1, 600)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occPop-1)){
    N[t+1] <- lambda * N[t]
  }

  # Observation model
  for (t in 1:n.occPop){
    C[t] ~ dnorm(N[t], tau)
  }
}
")


# 2. Definition of simulation parameters
# 2.1. Number of simulations
# nsim <- 1000    # ca 2.2 hours
nsim <- 5    # ~~~ for testing

# 2.2. Age specific survival probabilities (juv, adult)
phi <- c(0.3, 0.55)

# 2.3. Fecundity rate (females)
f <- 3.1

# 2.4. Initial population size per age class
Ni <- c(50, 50)

# 2.5. Number of years
T <- 10

# 2.6. Observation error for the population survey
sigma <- 10

# 2.7. Capture and recapture probabilities
cap <- c(0.4, 0.4)          # Capture prob. of nestlings and adults
recap <- 0.6                # Recapture probability

# 2.8. Probability to find a brood whose reproductive output is recorded
pprod <- 0.5

# 2.9. Define matrices to store the result
res1 <- array(NA, dim=c(45, 11, nsim))
res2 <- array(NA, dim=c(4, 11, nsim))
res3 <- array(NA, dim=c(2, 11, nsim))
res4 <- array(NA, dim=c(13, 11, nsim))
ss <- matrix(NA, nrow=nsim, ncol=4)       # For sample size

# 2.10. Settings for JAGS
# 2.10.1. Initial values
inits.ipm <- function(){list(mean.sj=runif(1, 0, 0.5))}
inits.cjs <- function(){list(mean.sj=runif(1, 0, 0.5))}
inits.pois <- function(){list(mean.f=runif(1, 1, 5))}
inits.ssm <- function(){list(sigma=runif(1, 2, 5))}

# 2.10.2. Parameters monitored
parameters.ipm <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma",
    "ann.growth.rate", "Ntot")
parameters.cjs <- c("mean.sj", "mean.sa", "mean.p")
parameters.pois <- c("mean.f")
parameters.ssm <- c("lambda", "N", "sigma")

# 2.10.3. MCMC settings
ni <- 20000; nb <- 10000; nc <- 3; nt <- 1; na <- 1000


# 3. Simulations
system.time(
for (s in 1:nsim){
  set.seed(s)
  # 3.1. Create 3 populations (independent)
  ind1 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
  ind2 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
  ind3 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)

  # 3.2. Create the population survey data
  count <- simCountNorm(N=ind1$totB, sigma=sigma)$count

  # 3.3. Create the capture histories and the corresponding m-arrays
  ch <- simCapHist(state=ind2$state, cap=cap, recap=recap, maxAge=2)
  marr <- marrayAge(ch$ch, ch$age)

  # 3.4. Create productivity data
  P <- simProd(reprod=ind3$reprod, pInclude=pprod)

  # Aggregate productivity data to make the model run faster
  sJ <- colSums(P$prod.agg)[1]
  nJ <- colSums(P$prod.agg)[2]

  # 3.5. Monitor sample size
  ss[s,1] <- mean(ind1$totA)         # Mean population size
  ss[s,2] <- table(ch$age)[1]        # Number of marked juveniles
  ss[s,3] <- table(ch$age)[2]        # Number of marked adults
  ss[s,4] <- nJ                      # Number of broods recorded

  # 3.6. Bundle data (4 sets)
  jags.data.ipm <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), sJ=sJ, nJ=nJ, C=count)
  jags.data.cjs <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occCJS=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]))
  jags.data.pois <- list(sJ=sJ, nJ=nJ)
  jags.data.ssm <- list(n.occPop=T, C=count)

  # 3.7. Call JAGS from R (jagsUI) to run the 4 models
  m1 <- try(jags(jags.data.ipm, inits.ipm, parameters.ipm, "model6.txt",
      n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE))
  if(!inherits(m1, "try-error"))
    res1[,,s] <- m1$summary

  m2 <- try(jags(jags.data.cjs, inits.cjs, parameters.cjs, "model7.txt",
      n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE))
  if(!inherits(m2, "try-error"))
    res2[,,s] <- m2$summary

  m3 <- try(jags(jags.data.pois, inits.pois, parameters.pois, "model8.txt",
      n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE))
  if(!inherits(m3, "try-error"))
    res3[,,s] <- m3$summary

  m4 <- try(jags(jags.data.ssm, inits.ssm, parameters.ssm, "model9.txt",
      n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE))
  if(!inherits(m4, "try-error"))
    res4[,,s] <- m4$summary

  print(s)
} ) #s

# 4. Save simulation results
save(res1, res2, res3, res4, phi, f, sigma, recap, ss, file="Data Fig 6.1.Rdata")

# 5. Produce figure 6.1
op <- par(cex=1.5)
boxplot(cbind(res1[1,2,]/res1[1,1,]*100, res2[1,2,]/res2[1,1,]*100,
    res1[2,2,]/res1[2,1,]*100, res2[2,2,]/res2[2,1,]*100, res1[4,2,]/res1[4,1,]*100,
    res3[1,2,]/res3[1,1,]*100, res1[34,2,]/res1[34,1,]*100, res4[1,2,]/res4[1,1,]*100),
    ylab="Coefficient of variation", ylim=c(0, 12), outline=FALSE,
    col=rep(c("red", "dodgerblue"), 4), border="black", axes=FALSE, boxwex=0.75, at=1:8)
axis(2, las=1)
axis(1, at = c(1.5, 3.5, 5.5, 7.5),
    labels = c(expression(italic('s')[italic(j)]),
    expression(italic('s')[italic(a)]), expression(italic('f')),
    expression(lambda)), tcl = -0.5, lwd = 1.5)
legend("topright", pch=rep(15,2), col=c("red", "dodgerblue"),
    legend=c("IPM", "Single data set"), bty="n")
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
