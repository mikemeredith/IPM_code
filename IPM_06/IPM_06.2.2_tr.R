# Schaub & Kéry (2022) Integrated Population Models
# Chapter 6 : Benefits of integrated population modeling
# ------------------------------------------------------
# Code from final MS.

# Run time testing 3 mins, full run 14 hrs

library(IPMbook) ; library(jagsUI)

# 6.2 Parameter estimates with increased precision
# ================================================

# 6.2.2 Where does the information come from?
# -------------------------------------------

# ~~~~ Code for the simulations: flow of information ~~~~

# 1. Analyzing models
# 1.1. Baseline IPM
cat(file="model10.txt", "
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



# 1.2. IPM with more CMR data
cat(file="model11.txt", "
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
  for (k in 1:ex){
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
      marr.j[t,1:n.occasions,k] ~ dmulti(pr.j[t,,k], rel.j[t,k])
      marr.a[t,1:n.occasions,k] ~ dmulti(pr.a[t,,k], rel.a[t,k])
    } #t
    # Define the cell probabilities of the m-arrays
    for (t in 1:(n.occasions-1)){
      # Main diagonal
      q[t,k] <- 1 - p[t]   # Probability of non-recapture
      pr.j[t,t,k] <- sj[t] * p[t]
      pr.a[t,t,k] <- sa[t] * p[t]
      # Above main diagonal
      for (j in (t+1):(n.occasions-1)){
         pr.j[t,j,k] <- sj[t] * prod(sa[(t+1):j]) * prod(q[t:(j-1),k]) * p[j]
         pr.a[t,j,k] <- prod(sa[t:j]) * prod(q[t:(j-1),k]) * p[j]
       } #j
       # Below main diagonal
       for (j in 1:(t-1)){
         pr.j[t,j,k] <- 0
         pr.a[t,j,k] <- 0
       } #j
     } #t
     # Last column: probability of non-recapture
     for (t in 1:(n.occasions-1)){
       pr.j[t,n.occasions,k] <- 1-sum(pr.j[t,1:(n.occasions-1),k])
       pr.a[t,n.occasions,k] <- 1-sum(pr.a[t,1:(n.occasions-1),k])
    } #t
  } #k

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


# 1.3. IPM with multiple productivity data
cat(file="model12.txt", "
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
  for (k in 1:ex){
    sJ[k] ~ dpois(nJ[k] * mean.f)
  }

  # Capture-recapture date (CJS model with multinomial likelihood)
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


# 1.4. IPM with more count data
cat(file="model13.txt", "
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
  for (k in 1:ex){
    for (t in 1:n.occasions){
      C[k,t] ~ dnorm(N[1,t] + N[2,t], tau)
    } #t
  } #k

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


# 2. Simulation parameters

# 2.1. Age specific survival probabilities (juv, adult)
phi <- c(0.3, 0.55)

# 2.2. Fecundity rate (females)
f <- 3.1

# 2.3. Initial population size per age class
Ni <- c(50, 50)

# 2.4. Number of years
T <- 10

# 2. 5. Observation error for the population survey
sigma <- 10

# 2.6. Capture and recapture probabilities
cap <- c(0.4, 0.4)          # Capture prob. of nestlings and adults
recap <- 0.6                # Recapture probability

# 2.7. Probability to find a brood whose reproductive output is recorded
pprod <- 0.5

# 2.8. Number of times data sets are replicated
ex <- 10

# 2.9. Number of simulations
# nsim <- 1300    # ca. 12.5 hours
nsim <- 4    # ~~~ testing

# 2.10. Define matrices to store the result
res1 <- res2 <- res3 <- res4 <- array(NA, dim=c(45, 11, nsim))

# 2.11. Settings for JAGS
# 2.11.1. Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5))}

# 2.11.2. Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma",
    "ann.growth.rate", "Ntot")

# 2.11.3. MCMC settings
ni <- 20000; nb <- 10000; nc <- 3; nt <- 1; na <- 2000


# 3. Start simulations
system.time(
for (s in 1:nsim){
  set.seed(s)
  # 3.1. Create 3 populations (to ensure independence)
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

  # 3.5. Bundle data
  # 3.5.1. Baseline: no replication
  jags.data.1 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), sJ=sJ, nJ=nJ, C=count)

  # 3.5.2. ex times more CMR data
  marr.j <- marr.a <- array(NA, dim=c(dim(marr)[1:2], ex))
  for (i in 1:ex){
    marr.j[,,i] <- marr[,,1]
    marr.a[,,i] <- marr[,,2]
  } #i
  jags.data.2 <- list(marr.j=marr.j, marr.a=marr.a, n.occasions=T,
      rel.j=apply(marr.j, c(1,3), sum), rel.a=apply(marr.a, c(1,3), sum),
      sJ=sJ, nJ=nJ, C=count, ex=ex)

  # 3.5.3. ex times more productivity data
  jags.data.3 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
      rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), sJ=rep(as.numeric(sJ), ex),
      nJ=rep(as.numeric(nJ), ex), C=count, ex=ex)

  # 3.5.4. ex times more count data
  C <- matrix(NA, nrow=ex, ncol=length(count))
  for (i in 1:ex){
    C[i,] <- count
  } #i
  jags.data.4 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), sJ=sJ, nJ=nJ, C=C, ex=ex)


  # 4. Run models in JAGS from R (jagsUI)
  # 4.1. Baseline
  m1 <- try(jags(jags.data.1, inits, parameters, "model10.txt",
      n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE))
  if(!inherits(m1, "try-error"))
    res1[,,s] <- m1$summary

  # 4.2. More CMR data
  m2 <- try(jags(jags.data.2, inits, parameters, "model11.txt",
      n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE))
  if(!inherits(m2, "try-error"))
    res2[,,s] <- m2$summary

  # 4.3. More productivity data
  m3 <- try(jags(jags.data.3, inits, parameters, "model12.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE))
  if(!inherits(m3, "try-error"))
    res3[,,s] <- m3$summary

  # 4.4. More count data
  m4 <- try(jags(jags.data.4, inits, parameters, "model13.txt",
      n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE))
  if(!inherits(m4, "try-error"))
    res4[,,s] <- m4$summary
  print(s)
} )#s

# 5. Save simulation results
save(res1, res2, res3, res4, phi, f, sigma, recap, ex, file="Data Fig 6.2.Rdata")

# 6. Produce Fig. 6.2
# Function to ensure that only converged estimates are included
r.incl <- function(res, param, r.crit=1.1){
  h <- which(res[param, 8, ] < r.crit)
  return(h)
}
co <- c("red", rep("dodgerblue", 3))
# op <- par(mfrow=c(2,2), cex=1.1, mar=c(3.5, 4, 1.5, 1))
op <- par(mfrow=c(2,2), mar=c(3.5, 4, 1.5, 1))
u <- which(as.numeric(table(c(r.incl(res1, 1), r.incl(res2, 1), r.incl(res3, 1),
    r.incl(res4, 1)))) == 4)[1:1000]
lab <- expression('CV ('*italic('s')[italic(j)]*')')
boxplot(cbind(res1[1,2,u] / res1[1,1,u] * 100, res2[1,2,u] / res2[1,1,u] * 100,
    res3[1,2,u] / res3[1,1,u] * 100, res4[1,2,u] / res4[1,1,u] * 100),
    ylab=lab, axes=FALSE, outline=FALSE, col=co)
axis(2, las=1)
axis(1, at=1:4, labels=NA)
mtext("Juvenile survival", at=2.5, line=0.5, adj=0.5, font=2)

u <- which(as.numeric(table(c(r.incl(res1, 2), r.incl(res2, 2), r.incl(res3, 2),
    r.incl(res4, 2)))) == 4)[1:1000]
lab <- expression('CV ('*italic('s')[italic(a)]*')')
boxplot(cbind(res1[2,2,u] / res1[2,1,u] * 100, res2[2,2,u] / res2[2,1,u] * 100,
    res3[2,2,u] / res3[2,1,u] * 100, res4[2,2,u] / res4[2,1,u] * 100), ylab=lab,
    axes=FALSE, outline=FALSE, col=co)
axis(2, las=1)
axis(1, at=1:4, labels=NA)
mtext("Adult survival", at=2.5, line=0.5, adj=0.5, font=2)

u <- which(as.numeric(table(c(r.incl(res1, 4), r.incl(res2, 4), r.incl(res3, 4),
    r.incl(res4, 4)))) == 4)[1:1000]
lab <- expression('CV ('*italic('f')*')')
boxplot(cbind(res1[4,2,u] / res1[4,1,u] * 100, res2[4,2,u] / res2[4,1,u] * 100,
    res3[4,2,u] / res3[4,1,u] * 100, res4[4,2,u] / res4[4,1,u] * 100), ylab=lab,
    axes=FALSE, outline=FALSE, col=co)
axis(2, las=1)
axis(1, at=1:4, labels=c("original", "10x CR", "10x Prod", "10x Count"))
mtext("Productivity", at=2.5, line=0.5, adj=0.5, font=2)

u <- which(as.numeric(table(c(r.incl(res1, 26), r.incl(res2, 26), r.incl(res3, 26),
    r.incl(res4, 26)))) == 4)[1:1000]
lab <- expression('CV ('*lambda*')')
boxplot(cbind(res1[26,2,u] / res1[26,1,u] * 100, res2[26,2,u] / res2[26,1,u] * 100,
    res3[26,2,u] / res3[26,1,u] * 100, res4[26,2,u] / res4[26,1,u] * 100), ylab=lab,
    axes=FALSE, outline=FALSE, col=co)
axis(2, las=1)
axis(1, at=1:4, labels=c("original", "10x CR", "10x Prod", "10x Count"))
mtext("Population growth rate", at=2.5, line=0.5, adj=0.5, font=2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
