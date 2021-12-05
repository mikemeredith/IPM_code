# Schaub & Kéry (2022) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------

# Run time approx. 2 mins

library(IPMbook) ; library(jagsUI)

# 4.5 Models for survival surveys
# ===============================

# 4.5.2 Multistate capture-recapture models
# -----------------------------------------

# 4.5.2.1 State-space formulation
# '''''''''''''''''''''''''''''''

# Choose constants in simulation
nmarked <- 50                           # Number of marked individuals each year and site
nyears <- 6                             # Number of years
phiA <- 0.8                             # Apparent survival probability when at site A
phiB <- 0.7                             # Apparent survival probability when at site B
psiAB <- 0.3                            # Probability to move from A to B
psiBA <- 0.5                            # Probability to move from B to A
pA <- 0.7                               # Recapture probability when at site A
pB <- 0.4                               # Recapture probability when at site B

# Determine occasion when an individual first captured and marked
f <- rep(rep(1:(nyears-1), each=nmarked), 2)
nind <- length(f)                       # Total number of marked individuals

# Construct the transition probability matrix (TPM), above called OMEGA
# Includes the dead state as state number 3
# Departure state in rows (time t), arrival state in columns (t+1)
TPM <- matrix(c(
    phiA * (1-psiAB), phiA * psiAB, 1-phiA,
    phiB * psiBA, phiB * (1-psiBA), 1-phiB,
    0, 0, 1 ), nrow=3, byrow=TRUE)

# Construct the observation probability matrix (OPM), above called THETA
# True state is in rows, observed state is in columns
# Includes nondetection as observation event number 3
OPM <- matrix(c(
    pA, 0, 1-pA,
    0, pB, 1-pB,
    0, 0, 1 ), nrow=3, byrow=TRUE)

# State or ecological process
# Simulate true system state
z <- array(NA, dim=c(nind, nyears))     # Empty alive/dead matrix

# Initial conditions: all individuals alive at f(i)
initial.state <- c(rep(1, nind/2), rep(2, nind/2))
for (i in 1:nind){
  z[i,f[i]] <- initial.state[i]
}

set.seed(2)                             # Initialize the RNGs in R
# Propagate alive/dead process forwards via transition rule (=TPM=OMEGA)
for (i in 1:nind){
  for (t in (f[i]+1):nyears){
    departure.state <- z[i,t-1]
    arrival.state <- which(rmultinom(1,1, TPM[departure.state,])==1)
    z[i,t] <- arrival.state
  } #t
} #i
z                                       # Not shown, but useful if you do look at this

# Observation process: simulate observations using observation matrix OPM (=THETA)
y <- array(3, dim=c(nind, nyears))
for (i in 1:nind){
  y[i,f[i]] <- z[i,f[i]]
  for (t in (f[i]+1):nyears){
    true.state <- z[i,t-1]
    observed.state <- which(rmultinom(1,1, OPM[true.state,])==1)
    y[i,t] <- observed.state
  } #t
} #i

head(y)
#      [,1] [,2] [,3] [,4] [,5] [,6]
# [1,]    1    1    1    3    3    3
# [2,]    1    1    3    3    3    3
# [3,]    1    1    1    3    3    3
# [4,]    1    3    1    3    3    3
# [5,]    1    1    2    2    1    1
# [6,]    1    3    3    3    3    3

# Data bundle
jags.data <- list(y=y, f=f, nind=nind, nyears=ncol(y))
str(jags.data)
# List of 4
# $ y     : num [1:500, 1:6] 1 1 1 1 1 1 1 1 1 1 ...
# $ f     : int [1:500] 1 1 1 1 1 1 1 1 1 1 ...
# $ nind  : int 500
# $ nyears: num 6

# Write JAGS model file
cat(file="model18.txt", "
model {
  # Priors and linear models
  for (t in 1:(nyears-1)){
    phiA[t] <- mean.phi[1]
    phiB[t] <- mean.phi[2]
    psiAB[t] <- mean.psi[1]
    psiBA[t] <- mean.psi[2]
    pA[t] <- mean.p[1]
    pB[t] <- mean.p[2]
  }
  for (u in 1:2){
    mean.phi[u] ~ dunif(0, 1)           # Priors for mean state-specific survival
    mean.psi[u] ~ dunif(0, 1)           # Priors for mean transitions
    mean.p[u] ~ dunif(0, 1)             # Priors for mean state-specific recapture
  }

  # Define state-transition and observation matrices
  for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    # (Transition probability matrix, or OMEGA)
    for (t in f[i]:(nyears-1)){
      TPM[1,i,t,1] <- phiA[t] * (1-psiAB[t])
      TPM[1,i,t,2] <- phiA[t] * psiAB[t]
      TPM[1,i,t,3] <- 1-phiA[t]
      TPM[2,i,t,1] <- phiB[t] * psiBA[t]
      TPM[2,i,t,2] <- phiB[t] * (1-psiBA[t])
      TPM[2,i,t,3] <- 1-phiB[t]
      TPM[3,i,t,1] <- 0
      TPM[3,i,t,2] <- 0
      TPM[3,i,t,3] <- 1

      # Define probabilities of observation O(t) given S(t)
      # (Observation probability matrix, or THETA)
      OPM[1,i,t,1] <- pA[t]
      OPM[1,i,t,2] <- 0
      OPM[1,i,t,3] <- 1-pA[t]
      OPM[2,i,t,1] <- 0
      OPM[2,i,t,2] <- pB[t]
      OPM[2,i,t,3] <- 1-pB[t]
      OPM[3,i,t,1] <- 0
      OPM[3,i,t,2] <- 0
      OPM[3,i,t,3] <- 1
    } #t
  } #i

  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):nyears){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(TPM[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(OPM[z[i,t], i, t-1,])
    } #t
  } #i
}
")

# Function to create a matrix of initial values for latent state z
zInitMS <- function(ch, f){
  states <- max(ch, na.rm=TRUE)
  known.states <- 1:(states-1)
  v <- which(ch==states)
  ch[v] <- sample(known.states, length(v), replace=TRUE)
  for (i in 1:nrow(ch)) ch[i,1:f[i]] <- NA
  return(ch)
}

# Initial values
inits <- function(){list(z=zInitMS(y, f))}

# Parameters monitored
parameters <- c("mean.phi", "mean.psi", "mean.p")

# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000

# Call JAGS from R (ART 3 min), check convergence and summarize posteriors
out21 <- jags(jags.data, inits, parameters, "model18.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out21) # Not shown
print(out21)
#               mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# mean.phi[1]  0.864  0.021    0.821    0.864    0.903    FALSE 1 1.002   957
# mean.phi[2]  0.813  0.028    0.756    0.814    0.866    FALSE 1 1.001  6000
# mean.psi[1]  0.120  0.019    0.086    0.119    0.159    FALSE 1 1.001  1503
# mean.psi[2]  0.230  0.026    0.183    0.229    0.283    FALSE 1 1.001  6000
# mean.p[1]    0.694  0.032    0.630    0.695    0.754    FALSE 1 1.004   885
# mean.p[2]    0.463  0.042    0.383    0.462    0.549    FALSE 1 1.001  2843


# 4.5.2.2 Multinomial formulation
# '''''''''''''''''''''''''''''''

y[y==3] <- 0

y[45:55,]
#       [,1] [,2] [,3] [,4] [,5] [,6]
#  [1,]    1    0    1    0    0    0
#  [2,]    1    1    1    0    0    0
#  [3,]    1    1    0    2    0    0
#  [4,]    1    0    1    1    0    0
#  [5,]    1    1    0    1    0    0
#  [6,]    1    1    1    2    1    1
#  [7,]    0    1    0    0    1    0
#  [8,]    0    1    0    0    0    0
#  [9,]    0    1    1    2    0    0
# [10,]    0    1    1    1    0    0
# [11,]    0    1    0    0    1    1

marr <- marray(y)
marr
#         recaptured
# released Y2.S1 Y2.S2 Y3.S1 Y3.S2 Y4.S1 Y4.S2 Y5.S1 Y5.S2 Y6.S1 Y6.S2 never
#    Y1.S1    34     0     8     1     1     1     0     0     1     0     4
#    Y1.S2     0    17     8     2     3     2     1     1     0     0    16
#    Y2.S1     0     0    51     3     9     2     3     1     0     1    14
#    Y2.S2     0     0     5    24     8     8     1     1     1     1    18
#    Y3.S1     0     0     0     0    56    10    13     5     3     1    34
#    Y3.S2     0     0     0     0     5    32     6     3     4     2    28
#    Y4.S1     0     0     0     0     0     0    72     3    14     3    40
#    Y4.S2     0     0     0     0     0     0    14    33    15     5    38
#    Y5.S1     0     0     0     0     0     0     0     0    85     9    66
#    Y5.S2     0     0     0     0     0     0     0     0    10    22    65

# Bundle data
# Calculate the number of states
ns <- length(unique(as.numeric(y))) - 1
jags.data <- list(marr=marr, nyears=ncol(y), rel=rowSums(marr), ns=ns, zero=matrix(0, ns, ns),
    ones=diag(ns))
str(jags.data)
# List of 6
# $ marr  : num [1:10, 1:11] 34 0 0 0 0 0 0 0 0 0 ...
# $ nyears: int 6
# $ rel   : num [1:10] 50 50 84 67 122 80 132 105 160 97
# $ ns    : num 2
# $ zero  : num [1:2, 1:2] 0 0 0 0
# $ ones  : num [1:2, 1:2] 1 0 0 1

# Write JAGS model file
cat(file="model19.txt", "
model {
  # Priors and linear models
  for (t in 1:(nyears-1)){
    phiA[t] <- mean.phi[1]
    phiB[t] <- mean.phi[2]
    psiAB[t] <- mean.psi[1]
    psiBA[t] <- mean.psi[2]
    pA[t] <- mean.p[1]
    pB[t] <- mean.p[2]
  }

  for (u in 1:2){
    mean.phi[u] ~ dunif(0, 1)           # Priors for mean state-specific survival
    mean.psi[u] ~ dunif(0, 1)           # Priors for mean transitions
    mean.p[u] ~ dunif(0, 1)             # Priors for mean state-specific recapture
  }

  # Define state-transition and re-encounter probabilities
  for (t in 1:(nyears-1)){
    psi[1,t,1] <- phiA[t] * (1-psiAB[t])
    psi[1,t,2] <- phiA[t] * psiAB[t]
    psi[2,t,1] <- phiB[t] * psiBA[t]
    psi[2,t,2] <- phiB[t] * (1-psiBA[t])

    po[1,t] <- pA[t]
    po[2,t] <- pB[t]
  }

  # From here onwards, no changes needed regardless of which model is fitted
  # Calculate probability of non-encounter (dq) and reshape the array for the encounter
  # probabilities
  for (t in 1:(nyears-1)){
    for (s in 1:ns){
      dp[s,t,s] <- po[s,t]
      dq[s,t,s] <- 1-po[s,t]
    } #s
    for (s in 1:(ns-1)){
      for (m in (s+1):ns){
        dp[s,t,m] <- 0
        dq[s,t,m] <- 0
      } #s
    } #m
    for (s in 2:ns){
      for (m in 1:(s-1)){
        dp[s,t,m] <- 0
        dq[s,t,m] <- 0
      } #s
    } #m
  } #t

  # Define the multinomial likelihood
  for (t in 1:((nyears-1)*ns)){
    marr[t,1:(nyears *ns-(ns-1))] ~ dmulti(pi[t,], rel[t])
  }

  # Define cell probabilities of the multistate m-array
  # Matrix U: product of probabilities of state-transition and non-encounter (needed because
  # there is no product function for matrix multiplication in JAGS)
  for (t in 1:(nyears-2)){
    U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- ones
    for (j in (t+1):(nyears-1)){
      U[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(j-2)*ns+(1:ns)] %*% psi[,t,] %*%
          dq[,t,]
    } #j
  } #t
  U[(nyears-2)*ns+(1:ns), (nyears-2)*ns+(1:ns)] <- ones

  for (t in 1:(nyears-2)){
    # Diagonal
    pi[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% psi[,t,] %*%
        dp[,t,]
    # Above main diagonal
    for (j in (t+1):(nyears-1)){
      pi[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] %*% psi[,j,] %*%
          dp[,j,]
    } #j
  } #t
  pi[(nyears-2)*ns+(1:ns),(nyears-2)*ns+(1:ns)] <- psi[,nyears-1,] %*% dp[,nyears-1,]

  # Below main diagonal
  for (t in 2:(nyears-1)){
    for (j in 1:(t-1)){
      pi[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
    } #j
  } #t

  # Last column: probability of non-recapture
  for (t in 1:((nyears-1)*ns)){
    pi[t,(nyears*ns-(ns-1))] <- 1-sum(pi[t,1:((nyears-1)*ns)])
  } #t
}
")

# Initial values
inits <- function(){list(mean.phi=runif(2, 0, 1))}

# Parameters monitored
parameters <- c("mean.phi", "mean.psi", "mean.p")

# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out22 <- jags(jags.data, inits, parameters, "model19.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
# par(mfrow=c(2, 2)); traceplot(out22) # Not shown
traceplot(out22) # Not shown
print(out22)
#                mean    sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# mean.phi[1]   0.863 0.021   0.822   0.863   0.902    FALSE 1 1.000  6000
# mean.phi[2]   0.813 0.029   0.756   0.813   0.867    FALSE 1 1.000  6000
# mean.psi[1]   0.119 0.020   0.085   0.118   0.161    FALSE 1 1.000  6000
# mean.psi[2]   0.230 0.027   0.180   0.229   0.285    FALSE 1 1.002  2185
# mean.p[1]     0.694 0.032   0.630   0.695   0.755    FALSE 1 1.000  6000
# mean.p[2]     0.466 0.043   0.385   0.465   0.553    FALSE 1 1.000  6000
