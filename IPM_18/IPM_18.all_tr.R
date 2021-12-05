# Schaub & Kéry (2022) Integrated Population Models
# Chapter 18 : Cormorant
# ----------------------

# Run time for test script 8 mins, full run 3 hrs

# 18.4 Component data likelihoods
# ===============================

library(IPMbook); library(jagsUI)
data(cormorant)
str(cormorant)
# List of 2
# $ count: int [1:3, 1:14] 5048 1982 804 4321 1860 1350 4634 2170 1848 4318 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:3] "V" "M" "S"
# .. ..$ : chr [1:14] "1991" "1992" "1993" "1994" ...
# $ ms.ch: int [1:12659, 1:14] 1 1 1 0 0 0 0 0 0 0 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:14] "1991" "1992" "1993" "1994" ...

marr <- marray(cormorant$ms.ch, unobs=3)


# 18.4.1 Population count data
# ----------------------------

phi <- c(0.4, 0.8)
kappa <- 0.25
rho <- 1.5
A <- matrix(c(
    phi[2] * (1-kappa), phi[1] * rho,
    phi[2] * kappa, phi[2]), byrow=TRUE, ncol=2)
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])
matrix(revec / sum(revec))                                # Standardized right eigenvector
#           [,1]
# [1,] 0.5657415
# [2,] 0.4342585


# 18.4.2 Multistate capture-recapture data (no code)


# 18.5 The integrated population model
# ====================================

# Bundle data and produce data overview
jags.data <- list(marr=marr, n.years=ncol(cormorant$ms.ch), rel=rowSums(marr), ns=9,
    zero=matrix(0, ncol=9, nrow=9), ones=diag(9), C=cormorant$count)
str(jags.data)
# List of 7
# $ marr   : num [1:117, 1:118] 0 0 0 0 0 0 0 0 0 0 ...
# $ n.years: int 14
# $ rel    : num [1:117] 379 0 0 532 39 0 0 0 0 375 ...
# $ ns     : num 9
# $ zero   : num [1:9, 1:9] 0 0 0 0 0 0 0 0 0 0 ...
# $ ones   : num [1:9, 1:9] 1 0 0 0 0 0 0 0 0 0 ...
# $ C      : int [1:3, 1:14] 5048 1982 804 4321 1860 1350 4634 2170 ...

# Write JAGS model file
cat(file = "model1.txt", "
model {
  # -------------------------------------------------
  # Stages:
  # N: not-yet recruited individuals
  # B: breeders
  # Parameters:
  # phi[age, site, time]: survival probability
  # eta[departure site, arrival site, time]: natal dispersal
  # nu[departure site, arrival site, time]: breeding dispersal
  # kappa[site, time]: recruitment probability
  # rho[site, time]: productivity
  # p[site, time]: recapture probability
  # -------------------------------------------------

  # Priors and linear models
  # Productivity
  for (t in 1:n.years){
    for (s in 1:3){
      rho[s,t] <- mean.rho[s]
    } #s
  } #t
  mean.rho[1] ~ dunif(0, 4)
  mean.rho[2] ~ dunif(0, 4)
  mean.rho[3] ~ dunif(0, 4)

  # Parameters of multistate model
  for (t in 1:(n.years-1)){
    phi[1,1,t] ~ dunif(0, 1)
    logit.phi[1,1,t] <- logit(phi[1,1,t])
    logit.phi[2,1,t] <- logit.phi[1,1,t] + beta.phi[2,1]
    phi[2,1,t] <- ilogit(logit.phi[2,1,t])
    logit.kappa[1,t] <- mu.kappa[1] + mu.kappa[2] * t
    kappa[1,t] <- ilogit(logit.kappa[1,t])
    for (s in 2:3){
      logit.phi[1,s,t] <- logit.phi[1,1,t] + beta.phi[1,s]
      phi[1,s,t] <- ilogit(logit.phi[1,s,t])
      logit.phi[2,s,t] <- logit.phi[1,1,t] + beta.phi[2,s]
      phi[2,s,t] <- ilogit(logit.phi[2,s,t])
      logit.kappa[s,t] <- logit.kappa[1,t] + mu.kappa[s+1]
      kappa[s,t] <- ilogit(logit.kappa[s,t])
    } #s
    for (s in 1:3){
      p[s,t] ~ dunif(0, 1)
      eta[1,s,t] <- mean.eta[1,s]
      eta[2,s,t] <- mean.eta[2,s]
      eta[3,s,t] <- mean.eta[3,s]
      nu[1,s,t] <- mean.nu[1,s]
      nu[2,s,t] <- mean.nu[2,s]
      nu[3,s,t] <- mean.nu[3,s]
    } #s
  } #t

  # Multinomial logit link to impose the constraint that natal dispersal does
  # only depend on the site of departure
  logit.mean.eta[1,2] ~ dnorm(0, 0.01)
  logit.mean.eta[1,3] <- logit.mean.eta[1,2]
  mean.eta[1,2] <- exp(logit.mean.eta[1,2]) / (1 + exp(logit.mean.eta[1,2]) +
      exp(logit.mean.eta[1,3]))
  mean.eta[1,3] <- exp(logit.mean.eta[1,3]) / (1 + exp(logit.mean.eta[1,2]) +
      exp(logit.mean.eta[1,3]))
  mean.eta[1,1] <- 1 - mean.eta[1,2] - mean.eta[1,3]
  logit.mean.eta[2,1] ~ dnorm(0, 0.01)
  logit.mean.eta[2,3] <- logit.mean.eta[2,1]
  mean.eta[2,1] <- exp(logit.mean.eta[2,1]) / (1 + exp(logit.mean.eta[2,1]) +
      exp(logit.mean.eta[2,3]))
  mean.eta[2,3] <- exp(logit.mean.eta[2,3]) / (1 + exp(logit.mean.eta[2,1]) +
      exp(logit.mean.eta[2,3]))
  mean.eta[2,2] <- 1 - mean.eta[2,1] - mean.eta[2,3]
  logit.mean.eta[3,1] ~ dnorm(0, 0.01)
  logit.mean.eta[3,2] <- logit.mean.eta[3,1]
  mean.eta[3,1] <- exp(logit.mean.eta[3,1]) / (1 + exp(logit.mean.eta[3,1]) +
      exp(logit.mean.eta[3,2]))
  mean.eta[3,2] <- exp(logit.mean.eta[3,2]) / (1 + exp(logit.mean.eta[3,1]) +
      exp(logit.mean.eta[3,2]))
  mean.eta[3,3] <- 1 - mean.eta[3,1] - mean.eta[3,2]

  # Multinomial logit link to impose the constraint that breeding dispersal does only depend on the
  # site of arrival
  logit.mean.nu[1,2] ~ dnorm(0, 0.01)
  logit.mean.nu[1,3] ~ dnorm(0, 0.01)
  logit.mean.nu[2,1] ~ dnorm(0, 0.01)
  mean.nu[1,2] <- exp(logit.mean.nu[1,2]) / (1 + exp(logit.mean.nu[1,2]) +
      exp(logit.mean.nu[1,3]) + exp(logit.mean.nu[2,1]))
  mean.nu[1,3] <- exp(logit.mean.nu[1,3]) / (1 + exp(logit.mean.nu[1,2]) +
      exp(logit.mean.nu[1,3]) + exp(logit.mean.nu[2,1]))
  mean.nu[2,1] <- exp(logit.mean.nu[2,1]) / (1 + exp(logit.mean.nu[1,2]) +
      exp(logit.mean.nu[1,3]) + exp(logit.mean.nu[2,1]))
  mean.nu[1,1] <- 1 - mean.nu[1,2] - mean.nu[1,3]
  mean.nu[2,2] <- 1 - mean.nu[2,1] - mean.nu[1,3]
  mean.nu[2,3] <- mean.nu[1,3]
  mean.nu[3,1] <- mean.nu[2,1]
  mean.nu[3,2] <- mean.nu[1,2]
  mean.nu[3,3] <- 1 - mean.nu[3,1] - mean.nu[3,2]

  for (i in 1:4){
    mu.kappa[i] ~ dnorm(0, 0.01)
  }

  beta.phi[2,1] ~ dnorm(0, 0.01)
  beta.phi[1,2] ~ dnorm(0, 0.01)
  beta.phi[2,2] ~ dnorm(0, 0.01)
  beta.phi[1,3] ~ dnorm(0, 0.01)
  beta.phi[2,3] ~ dnorm(0, 0.01)

  # Residual (observation) error
  for (s in 1:3){
    sigma[s] ~ dunif(0.05, 1000)
    tau[s] <- pow(sigma[s], -2)
  }

  # Population count data (state-space model)
  # Models for the initial population size: uniform priors
  N[1,1] ~ dunif(4270, 7940)
  B[1,1] ~ dunif(3500, 6500)
  N[2,1] ~ dunif(1710, 3180)
  B[2,1] ~ dunif(1400, 2600)
  N[3,1] ~ dunif(680, 1270)
  B[3,1] ~ dunif(560, 1040)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.years-1)){
    N[1,t+1] <- B[1,t] * rho[1,t] * phi[1,1,t] * eta[1,1,t] + B[2,t] * rho[2,t] * phi[1,2,t] *
        eta[2,1,t] + B[3,t] * rho[3,t] * phi[1,3,t] * eta[3,1,t] + N[1,t] * phi[2,1,t] *
        (1 - kappa[1,t])
    N[2,t+1] <- B[1,t] * rho[1,t] * phi[1,1,t] * eta[1,2,t] + B[2,t] * rho[2,t] * phi[1,2,t] *
        eta[2,2,t] + B[3,t] * rho[3,t] * phi[1,3,t] * eta[3,2,t] + N[2,t] * phi[2,2,t] *
        (1 - kappa[2,t])
    N[3,t+1] <- B[1,t] * rho[1,t] * phi[1,1,t] * eta[1,3,t] + B[2,t] * rho[2,t] * phi[1,2,t] *
        eta[2,3,t] + B[3,t] * rho[3,t] * phi[1,3,t] * eta[3,3,t] + N[3,t] * phi[2,3,t] *
        (1 - kappa[3,t])
    B[1,t+1] <- B[1,t] * phi[2,1,t] * nu[1,1,t] + B[2,t] * phi[2,2,t] * nu[2,1,t] + B[3,t] *
        phi[2,3,t] * nu[3,1,t] + N[1,t] * phi[2,1,t] * kappa[1,t]
    B[2,t+1] <- B[1,t] * phi[2,1,t] * nu[1,2,t] + B[2,t] * phi[2,2,t] * nu[2,2,t] + B[3,t] *
        phi[2,3,t] * nu[3,2,t] + N[2,t] * phi[2,2,t] * kappa[2,t]
    B[3,t+1] <- B[1,t] * phi[2,1,t] * nu[1,3,t] + B[2,t] * phi[2,2,t] * nu[2,3,t] + B[3,t] *
        phi[2,3,t] * nu[3,3,t] + N[3,t] * phi[2,3,t] * kappa[3,t]
  }

  # Observation model
  for (t in 1:n.years){
    for (s in 1:3){
      C[s,t] ~ dnorm(B[s,t], tau[s])
    } #s
  } #t

  # Multistate capture-recapture data (with multinomial likelihood)
  # Define state-transition and re-encounter probabilities
  for (t in 1:(n.years-1)){
    psi[1,t,1] <- 0
    psi[1,t,2] <- 0
    psi[1,t,3] <- 0
    psi[1,t,4] <- 0
    psi[1,t,5] <- 0
    psi[1,t,6] <- 0
    psi[1,t,7] <- phi[1,1,t] * eta[1,1,t]
    psi[1,t,8] <- phi[1,1,t] * eta[1,2,t]
    psi[1,t,9] <- phi[1,1,t] * eta[1,3,t]
    psi[2,t,1] <- 0
    psi[2,t,2] <- 0
    psi[2,t,3] <- 0
    psi[2,t,4] <- 0
    psi[2,t,5] <- 0
    psi[2,t,6] <- 0
    psi[2,t,7] <- phi[1,2,t] * eta[2,1,t]
    psi[2,t,8] <- phi[1,2,t] * eta[2,2,t]
    psi[2,t,9] <- phi[1,2,t] * eta[2,3,t]
    psi[3,t,1] <- 0
    psi[3,t,2] <- 0
    psi[3,t,3] <- 0
    psi[3,t,4] <- 0
    psi[3,t,5] <- 0
    psi[3,t,6] <- 0
    psi[3,t,7] <- phi[1,3,t] * eta[3,1,t]
    psi[3,t,8] <- phi[1,3,t] * eta[3,2,t]
    psi[3,t,9] <- phi[1,3,t] * eta[3,3,t]
    psi[4,t,1] <- 0
    psi[4,t,2] <- 0
    psi[4,t,3] <- 0
    psi[4,t,4] <- phi[2,1,t] * nu[1,1,t]
    psi[4,t,5] <- phi[2,1,t] * nu[1,2,t]
    psi[4,t,6] <- phi[2,1,t] * nu[1,3,t]
    psi[4,t,7] <- 0
    psi[4,t,8] <- 0
    psi[4,t,9] <- 0
    psi[5,t,1] <- 0
    psi[5,t,2] <- 0
    psi[5,t,3] <- 0
    psi[5,t,4] <- phi[2,2,t] * nu[2,1,t]
    psi[5,t,5] <- phi[2,2,t] * nu[2,2,t]
    psi[5,t,6] <- phi[2,2,t] * nu[2,3,t]
    psi[5,t,7] <- 0
    psi[5,t,8] <- 0
    psi[5,t,9] <- 0
    psi[6,t,1] <- 0
    psi[6,t,2] <- 0
    psi[6,t,3] <- 0
    psi[6,t,4] <- phi[2,3,t] * nu[3,1,t]
    psi[6,t,5] <- phi[2,3,t] * nu[3,2,t]
    psi[6,t,6] <- phi[2,3,t] * nu[3,3,t]
    psi[6,t,7] <- 0
    psi[6,t,8] <- 0
    psi[6,t,9] <- 0
    psi[7,t,1] <- 0
    psi[7,t,2] <- 0
    psi[7,t,3] <- 0
    psi[7,t,4] <- phi[2,1,t] * kappa[1,t]
    psi[7,t,5] <- 0
    psi[7,t,6] <- 0
    psi[7,t,7] <- phi[2,1,t] * (1-kappa[1,t])
    psi[7,t,8] <- 0
    psi[7,t,9] <- 0
    psi[8,t,1] <- 0
    psi[8,t,2] <- 0
    psi[8,t,3] <- 0
    psi[8,t,4] <- 0
    psi[8,t,5] <- phi[2,2,t] * kappa[2,t]
    psi[8,t,6] <- 0
    psi[8,t,7] <- 0
    psi[8,t,8] <- phi[2,2,t] * (1-kappa[2,t])
    psi[8,t,9] <- 0
    psi[9,t,1] <- 0
    psi[9,t,2] <- 0
    psi[9,t,3] <- 0
    psi[9,t,4] <- 0
    psi[9,t,5] <- 0
    psi[9,t,6] <- phi[2,3,t] * kappa[3,t]
    psi[9,t,7] <- 0
    psi[9,t,8] <- 0
    psi[9,t,9] <- phi[2,3,t] * (1-kappa[3,t])

    po[1,t] <- 0
    po[2,t] <- 0
    po[3,t] <- 0
    po[4,t] <- p[1,t]
    po[5,t] <- p[2,t]
    po[6,t] <- p[3,t]
    po[7,t] <- 0
    po[8,t] <- 0
    po[9,t] <- 0

    # Calculate probability of non-encounter (dq) and reshape the array for the encounter
    # probabilities
    for (s in 1:ns){
      dp[s,t,s] <- po[s,t]
      dq[s,t,s] <- 1-po[s,t]
    } #s
    for (s in 1:(ns-1)){
      for (m in (s+1):ns){
        dp[s,t,m] <- 0
        dq[s,t,m] <- 0
      } #m
    } #s
    for (s in 2:ns){
      for (m in 1:(s-1)){
        dp[s,t,m] <- 0
        dq[s,t,m] <- 0
      } #m
    } #s
  } #t

  # Define the multinomial likelihood
  for (t in 1:((n.years-1)*ns)){
    marr[t,1:(n.years*ns-(ns-1))] ~ dmulti(pr[t,], rel[t])
  }

  # Define the cell probabilities of the multistate m-array
  for (t in 1:(n.years-2)){
    U[(t-1)*ns+(1:ns), (t-1)*ns+(1:ns)] <- ones
    for (j in (t+1):(n.years-1)){
      U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*%
          psi[,t,] %*% dq[,t,]
    } #j
  } #t
  U[(n.years-2)*ns+(1:ns), (n.years-2)*ns+(1:ns)] <- ones

  # Diagonal
  for (t in 1:(n.years-2)){
    pr[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)]
        %*% psi[,t,] %*% dp[,t,]
    # Above main diagonal
    for (j in (t+1):(n.years-1)){
      pr[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] %*%
          psi[,j,] %*% dp[,j,]
    } #j
  } #t
  pr[(n.years-2)*ns+(1:ns), (n.years-2)*ns+(1:ns)] <- psi[,n.years-1,] %*% dp[,n.years-1,]

  # Below main diagonal
  for (t in 2:(n.years-1)){
    for (j in 1:(t-1)){
      pr[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
    } #j
  } #t

  # Last column: probability of non-recapture
  for (t in 1:((n.years-1)*ns)){
    pr[t,(n.years*ns-(ns-1))] <- 1-sum(pr[t,1:((n.years-1)*ns)])
  } #t
}
")

# Initial values
inits <- function(cc=cormorant$count){
  B <- array(NA, dim=c(3, 14))
  B[1,1] <- rpois(1, cc[1,1])
  B[2,1] <- rpois(1, cc[2,1])
  B[3,1] <- rpois(1, cc[3,1])
  N <- B * 0.55/0.45
  return(list(B=B, N=N))
}

# Parameters monitored
parameters <- c("beta.phi", "phi", "mu.kappa", "kappa", "p", "mean.eta", "mean.nu", "mean.rho", "sigma",
    "N", "B")

# MCMC settings
# ni <- 150000; nb <- 50000; nc <- 3; nt <- 100; na <- 3000
ni <- 3000; nb <- 1000; nt <- 1; nc <- 3; na <- 3000  # ~~~ for testing

# Call JAGS from R (ART 143 min) and check convergence
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1)


# 18.6. Results
# =============

print(out1, 3)
#                    mean       sd     2.5%       50%     97.5% overlap0     f  Rhat n.eff
# mu.phi[2,1]       0.201    0.298   -0.410     0.212     0.753     TRUE 0.765 1.016   167
# mu.phi[1,2]      -1.654    0.333   -2.320    -1.652    -0.999    FALSE 1.000 1.011   223
# mu.phi[2,2]       0.244    0.291   -0.347     0.258     0.779     TRUE 0.802 1.016   171
# mu.phi[1,3]      -0.648    0.432   -1.470    -0.658     0.267     TRUE 0.936 1.002   834
# mu.phi[2,3]       0.452    0.300   -0.167     0.456     1.017     TRUE 0.929 1.019   138
# phi[1,1,1]        0.686    0.063    0.562     0.687     0.807    FALSE 1.000 1.017   140
# phi[2,1,1]        0.731    0.016    0.698     0.731     0.761    FALSE 1.000 1.001  3000
# ...[output truncated]...
# phi[1,3,13]       0.582    0.094    0.416     0.577     0.778    FALSE 1.000 1.000  2855
# phi[2,3,13]       0.806    0.040    0.737     0.802     0.895    FALSE 1.000 1.002  1077
# mu.kappa[1]      -0.774    0.189   -1.118    -0.778    -0.381    FALSE 1.000 1.006   325
# mu.kappa[2]      -0.160    0.019   -0.198    -0.160    -0.123    FALSE 1.000 1.000  3000
# mu.kappa[3]       2.135    0.301    1.525     2.141     2.706    FALSE 1.000 1.003   713
# mu.kappa[4]       1.336    0.249    0.846     1.330     1.834    FALSE 1.000 1.001  1486
# kappa[1,1]        0.283    0.036    0.221     0.281     0.363    FALSE 1.000 1.007   286
# kappa[2,1]        0.764    0.054    0.645     0.769     0.855    FALSE 1.000 1.001  3000
# ...[output truncated]...
# kappa[2,13]       0.329    0.056    0.228     0.328     0.443    FALSE 1.000 1.000  3000
# kappa[3,13]       0.181    0.033    0.123     0.178     0.252    FALSE 1.000 1.000  3000
# p[1,1]            0.648    0.021    0.607     0.648     0.688    FALSE 1.000 1.000  3000
# p[2,1]            0.036    0.028    0.001     0.029     0.105    FALSE 1.000 1.001  1883
# ...[output truncated]...
# p[2,13]           0.947    0.042    0.841     0.956     0.998    FALSE 1.000 1.001  2164
# p[3,13]           0.650    0.051    0.552     0.651     0.751    FALSE 1.000 1.003   637
# mean.eta[1,1]     0.978    0.005    0.966     0.979     0.987    FALSE 1.000 1.000  3000
# mean.eta[2,1]     0.086    0.018    0.054     0.085     0.124    FALSE 1.000 1.002   932
# mean.eta[3,1]     0.005    0.002    0.001     0.004     0.011    FALSE 1.000 1.001  2490
# mean.eta[1,2]     0.011    0.003    0.006     0.011     0.017    FALSE 1.000 1.000  3000
# mean.eta[2,2]     0.827    0.035    0.752     0.829     0.892    FALSE 1.000 1.002   932
# mean.eta[3,2]     0.005    0.002    0.001     0.004     0.011    FALSE 1.000 1.001  2490
# mean.eta[1,3]     0.011    0.003    0.006     0.011     0.017    FALSE 1.000 1.000  3000
# mean.eta[2,3]     0.086    0.018    0.054     0.085     0.124    FALSE 1.000 1.002   932
# mean.eta[3,3]     0.991    0.005    0.979     0.992     0.998    FALSE 1.000 1.001  2490
# mean.nu[1,1]      0.989    0.002    0.984     0.989     0.993    FALSE 1.000 1.001  1674
# mean.nu[2,1]      0.006    0.002    0.003     0.006     0.011    FALSE 1.000 1.002  1362
# mean.nu[3,1]      0.006    0.002    0.003     0.006     0.011    FALSE 1.000 1.002  1362
# mean.nu[1,2]      0.005    0.001    0.003     0.005     0.008    FALSE 1.000 1.001  3000
# mean.nu[2,2]      0.988    0.003    0.982     0.988     0.993    FALSE 1.000 1.002  1139
# mean.nu[3,2]      0.005    0.001    0.003     0.005     0.008    FALSE 1.000 1.001  3000
# mean.nu[1,3]      0.006    0.002    0.003     0.006     0.010    FALSE 1.000 1.001  1898
# mean.nu[2,3]      0.006    0.002    0.003     0.006     0.010    FALSE 1.000 1.001  1898
# mean.nu[3,3]      0.989    0.002    0.983     0.989     0.993    FALSE 1.000 1.001  1398
# mean.rho[1]       1.650    0.133    1.390     1.647     1.908    FALSE 1.000 1.001  2853
# mean.rho[2]       2.220    0.204    1.840     2.212     2.630    FALSE 1.000 1.001  1994
# mean.rho[3]       1.641    0.159    1.347     1.631     1.966    FALSE 1.000 1.001  1301
# sigma[1]        303.682   91.371  171.488   289.210   526.688    FALSE 1.000 1.002   839
# sigma[2]        265.069   71.951  163.375   251.908   441.373    FALSE 1.000 1.000  3000
# sigma[3]        320.671   82.020  202.394   306.157   524.126    FALSE 1.000 1.000  3000
# N[1,1]         4902.564  595.877 4288.845  4723.980  6474.122    FALSE 1.000 1.000  3000
# N[2,1]         1867.905  150.467 1713.960  1821.940  2269.922    FALSE 1.000 1.001  3000
# ...[output truncated]...
# N[2,14]        2237.194  516.186 1416.561  2178.520  3424.011    FALSE 1.000 1.001  3000
# N[3,14]        7498.992 1512.295 4967.528  7374.228 10861.533    FALSE 1.000 1.000  3000
# B[1,1]         4678.304  270.860 4089.630  4705.478  5148.747    FALSE 1.000 1.002  1119
# B[2,1]         1604.252  132.981 1409.873  1586.594  1901.241    FALSE 1.000 1.000  3000
# ...[output truncated]...
# B[2,14]        1931.980  153.149 1640.733  1928.051  2237.172    FALSE 1.000 1.000  3000
# B[3,14]        3564.089  211.356 3187.816  3552.110  4016.581    FALSE 1.000 1.000  2507


# ~~~~ Code for Fig. 18.3 ~~~~
library(RColorBrewer)
co1 <- brewer.pal(n = 8, name = 'Blues')[c(8,5)]
co2 <- brewer.pal(n = 8, name = 'Greens')[c(5,8)]

op <- par(mar=c(1, 4.5, 2, 1), "mfrow")
layout(matrix(1:6, 3, 2, byrow=TRUE), widths=c(1.6, 1.6), heights=c(1,1,1.15), TRUE)
years <- 1991:2004
ny <- length(years)-1

site <- 1
d <- 0.15
plot(y=out1$mean$B[site,], x=(1:length(years))-d, type="b", pch=16, ylim=c(0, 14000),
    axes=FALSE, ylab="Population size", xlab=NA, col=co1[1])
mtext("Vorsø", side=3, line=-0.25)
segments((1:length(years))-d, out1$q2.5$B[site,], (1:length(years))-d,
    out1$q97.5$B[site,], col=co1[1])
axis(1, at=1:length(years), tcl=-0.25, labels=NA)
axis(1, at=seq(1, length(years), by=2), tcl=-0.5, labels=NA)
axis(2, at=c(0, 2000, 4000, 6000, 8000, 10000, 12000), tcl=-0.25, labels=NA)
axis(2, las=1, at=c(0, 4000, 8000, 12000), tcl=-0.5, labels=c(0, 4000, 8000, 12000))
points(y=out1$mean$N[site,], x=(1:length(years))+d, type="b", pch=16, col=co1[2])
segments((1:length(years))+d, out1$q2.5$N[site,], (1:length(years))+d, out1$q97.5$N[site,],
    col=co1[2])
points(cormorant$count[site,], pch=3, col="red")

plot(y=out1$mean$phi[1,site,], x=(1:ny)+0.5-d, type="b", pch=16, ylim=c(0, 1), axes=FALSE,
    ylab="Survival probability", xlab=NA, xlim=c(1,length(years)), col=co2[1])
segments((1:ny)+0.5-d, out1$q2.5$phi[1,site,], (1:ny)+0.5-d,
    out1$q97.5$phi[1,site,], col=co2[1])
axis(1, at=1:length(years), tcl=-0.25, labels=NA)
axis(1, at=seq(1, length(years), by=2), tcl=-0.5, labels=NA)
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), tcl=-0.25, labels=NA)
axis(2, las=1, at=c(0, 0.4, 0.8), tcl=-0.5, labels=c("0", "0.4", "0.8"))
points(y=out1$mean$phi[2,site,], x=(1:ny)+0.5+d, type="b", pch=16, col=co2[2])
segments((1:ny)+0.5+d, out1$q2.5$phi[2,site,], (1:ny)+0.5+d,
    out1$q97.5$phi[2,site,], col=co2[2])
legend("bottomleft", pch=c(16,16), col=co2, legend=c("Juvenile", "Adult"), bty="n")

site <- 2
d <- 0.15
plot(y=out1$mean$B[site,], x=(1:length(years))-d, type="b", pch=16, ylim=c(0, 4000),
    axes=FALSE, ylab="Population size", xlab=NA, col=co1[1])
mtext("Mågeøerne", side=3, line=-0.25)
segments((1:length(years))-d, out1$q2.5$B[site,], (1:length(years))-d,
    out1$q97.5$B[site,], col=co1[1])
axis(1, at=1:length(years), tcl=-0.25, labels=NA)
axis(1, at=seq(1, length(years), by=2), tcl=-0.5, labels=NA)
axis(2, at=c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500), tcl=-0.25, labels=NA)
axis(2, las=1, at=c(0, 1000, 2000, 3000, 4000), tcl=-0.5, labels=c(0, 1000, 2000, 3000, 4000))
points(y=out1$mean$N[site,], x=(1:length(years))+d, type="b", pch=16, col=co1[2])
segments((1:length(years))+d, out1$q2.5$N[site,], (1:length(years))+d,
    out1$q97.5$N[site,], col=co1[2])
points(cormorant$count[site,], pch=3, col="red")

plot(y=out1$mean$phi[1,site,], x=(1:ny)+0.5-d, type="b", pch=16, ylim=c(0, 1), axes=FALSE,
    ylab="Survival probability", xlab=NA, xlim=c(1,length(years)), col=co2[1])
segments((1:ny)+0.5-d, out1$q2.5$phi[1,site,], (1:ny)+0.5-d,
    out1$q97.5$phi[1,site,], col=co2[1])
axis(1, at=1:length(years), tcl=-0.25, labels=NA)
axis(1, at=seq(1, length(years), by=2), tcl=-0.5, labels=NA)
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), tcl=-0.25, labels=NA)
axis(2, las=1, at=c(0, 0.4, 0.8), tcl=-0.5, labels=c("0", "0.4", "0.8"))
points(y=out1$mean$phi[2,site,], x=(1:ny)+0.5+d, type="b", pch=16, col=co2[2])
segments((1:ny)+0.5+d, out1$q2.5$phi[2,site,], (1:ny)+0.5+d,
    out1$q97.5$phi[2,site,], col=co2[2])

par(mar=c(3, 4.5, 2, 1))
site <- 3
d <- 0.15
plot(y=out1$mean$B[site,], x=(1:length(years))-d, type="b", pch=16, ylim=c(0, 12000),
    axes=FALSE, ylab="Population size", xlab=NA, col=co1[1])
segments((1:length(years))-d, out1$q2.5$B[site,], (1:length(years))-d,
    out1$q97.5$B[site,], col=co1[1])
mtext("Stavns Fjord", side=3, line=-0.25)
axis(1, at=1:length(years), tcl=-0.25, labels=NA)
axis(1, at=seq(1, length(years), by=2), tcl=-0.5, labels=seq(years[1], years[14], by=2))
axis(2, at=c(0, 2000, 4000, 6000, 8000, 10000, 12000), tcl=-0.25, labels=NA)
axis(2, las=1, at=c(0, 4000, 8000, 12000), tcl=-0.5, labels=c(0, 4000, 8000, 12000))
points(y=out1$mean$N[site,], x=(1:length(years))+d, type="b", pch=16, col=co1[2])
segments((1:length(years))+d, out1$q2.5$N[site,], (1:length(years))+d,
    out1$q97.5$N[site,], col=co1[2])
points(cormorant$count[site,], pch=3, col="red")
legend("topleft", pch=c(16, 16, 3), col=c(co1, "red"),
    legend=c("Breeders", "Non-breeders", "Counts"), bty="n")

plot(y=out1$mean$phi[1,site,], x=(1:ny)+0.5-d, type="b", pch=16, ylim=c(0, 1),
    axes=FALSE, ylab="Survival probability", xlab=NA, xlim=c(1,length(years)), col=co2[1])
segments((1:ny)+0.5-d, out1$q2.5$phi[1,site,], (1:ny)+0.5-d,
    out1$q97.5$phi[1,site,], col=co2[1])
axis(1, at=1:length(years), tcl=-0.25, labels=NA)
axis(1, at=seq(1, length(years), by=2), tcl=-0.5, labels=seq(years[1], years[14], by=2))
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), tcl=-0.25, labels=NA)
axis(2, las=1, at=c(0, 0.4, 0.8), tcl=-0.5, labels=c("0", "0.4", "0.8"))
points(y=out1$mean$phi[2,site,], x=(1:ny)+0.5+d, type="b", pch=16, col=co2[2])
segments((1:ny)+0.5+d, out1$q2.5$phi[2,site,], (1:ny)+0.5+d,
    out1$q97.5$phi[2,site,], col=co2[2])
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Code for Fig. 18.4 ~~~~
library(denstrip)
library(plotrix)
op <- par(mar=c(4.5, 6, 3, 1), "mfrow")
layout(matrix(1:3,3,1,byrow=TRUE), widths=1.6, heights=c(1,1,1), TRUE)
plot(0, ylim=c(0.8, 3.2), xlim=c(0.725,1), axes=FALSE, pch=NA,
    xlab="Natal site fidelity", ylab=NA)
denstrip(out1$sims.list$mean.eta[,1,1], at=3,
    ticks=c(out1$mean$mean.eta[1,1], out1$q2.5$mean.eta[1,1], out1$q97.5$mean.eta[1,1]),
    twd=c(4,1.5,1.5), tlen=c(2,2,2), width=1/5)
denstrip(out1$sims.list$mean.eta[,2,2], at=2,
    ticks=c(out1$mean$mean.eta[2,2], out1$q2.5$mean.eta[2,2], out1$q97.5$mean.eta[2,2]),
    twd=c(4,1.5,1.5), tlen=c(2,2,2), width=1/5)
denstrip(out1$sims.list$mean.eta[,3,3], at=1,
    ticks=c(out1$mean$mean.eta[3,3], out1$q2.5$mean.eta[3,3], out1$q97.5$mean.eta[3,3]),
    twd=c(4,1.5,1.5), tlen=c(2,2,2), width=1/5)
axis(1)
axis(2, las=1, at=c(3,2,1), labels=c("Vorsø", "Mågeøerne", "Stavns\n Fjord"))
corner.label('A', font=2, cex=1.25)

plot(0, ylim=c(0.8, 3.2), xlim=c(0.973,1), axes=FALSE, pch=NA,
    xlab="Breeding site fidelity", ylab=NA)
denstrip(out1$sims.list$mean.nu[,1,1], at=3,
    ticks=c(out1$mean$mean.nu[1,1], out1$q2.5$mean.nu[1,1], out1$q97.5$mean.nu[1,1]),
    twd=c(4,1.5,1.5), tlen=c(2,2,2), width=1/5)
denstrip(out1$sims.list$mean.nu[,2,2], at=2,
    ticks=c(out1$mean$mean.nu[2,2], out1$q2.5$mean.nu[2,2], out1$q97.5$mean.nu[2,2]),
    twd=c(4,1.5,1.5), tlen=c(2,2,2), width=1/5)
denstrip(out1$sims.list$mean.nu[,3,3], at=1,
    ticks=c(out1$mean$mean.nu[3,3], out1$q2.5$mean.nu[3,3], out1$q97.5$mean.nu[3,3]),
    twd=c(4,1.5,1.5), tlen=c(2,2,2), width=1/5)
axis(1, at=c(0.975, 0.98, 0.985, 0.99, 0.995, 1), labels=NA, tcl=-0.25)
axis(1, at=c(0.98, 0.99, 1), labels=c('0.98', '0.99', '1.00'))

axis(2, las=1, at=c(3,2,1), labels=c("Vorsø", "Mågeøerne", "Stavns\n Fjord"))
corner.label('B', font=2, cex=1.25)

plot(0, ylim=c(0.8, 3.2), xlim=c(1,3), axes=FALSE, pch=NA, xlab="Productivity", ylab=NA)
denstrip(out1$sims.list$mean.rho[,1], at=3,
    ticks=c(out1$mean$mean.rho[1], out1$q2.5$mean.rho[1], out1$q97.5$mean.rho[1]),
    twd=c(4,1.5,1.5), tlen=c(2,2,2), width=1/5)
denstrip(out1$sims.list$mean.rho[,2], at=2,
    ticks=c(out1$mean$mean.rho[2], out1$q2.5$mean.rho[2], out1$q97.5$mean.rho[2]),
    twd=c(4,1.5,1.5), tlen=c(2,2,2), width=1/5)
denstrip(out1$sims.list$mean.rho[,3], at=1,
    ticks=c(out1$mean$mean.rho[3], out1$q2.5$mean.rho[3], out1$q97.5$mean.rho[3]),
    twd=c(4,1.5,1.5), tlen=c(2,2,2), width=1/5)
axis(1)
axis(2, las=1, at=c(3,2,1), labels=c("Vorsø", "Mågeøerne", "Stavns\n Fjord"))
corner.label('C', font=2, cex=1.25)
par(op)

# ~~~~ Code for Fig. 18.5 ~~~~
library(scales)
co <- viridis_pal(option='E')(20)[c(1,10,17)]
years <- 1991:2004
ny <- length(years)-1

op <- par(mar=c(3.5, 4.5, 1, 1))
plot(y=out1$mean$kappa[1,], x=(1:ny)+0.5, type="b", pch=16, ylim=c(0, 1), axes=FALSE,
    ylab='Recruitment probability', xlab=NA, xlim=c(1,length(years)), col=co[1])
segments((1:ny)+0.5, out1$q2.5$kappa[1,], (1:ny)+0.5, out1$q97.5$kappa[1,], col=co[1])
axis(1, at=1:length(years), tcl=-0.25, labels=NA)
axis(1, at=seq(1, length(years), by=2), tcl=-0.5, labels=seq(years[1], years[14], by=2))
axis(2, at=c(0.1, 0.3, 0.5, 0.7, 0.9), tcl=-0.25, labels=NA)
axis(2, las=1, at=c(0, 0.2, 0.4, 0.6, 0.8), tcl=-0.5, labels=c(0, 0.2, 0.4, 0.6, 0.8))
points(y=out1$mean$kappa[2,], x=(1:ny)+0.5+d, type="b", pch=16, col=co[2])
segments((1:ny)+0.5+d, out1$q2.5$kappa[2,], (1:ny)+0.5+d, out1$q97.5$kappa[2,], col=co[2])
points(y=out1$mean$kappa[3,], x=(1:ny)+0.5-d, type="b", pch=16, col=co[3])
segments((1:ny)+0.5-d, out1$q2.5$kappa[3,], (1:ny)+0.5-d, out1$q97.5$kappa[3,], col=co[3])
legend('topright', pch=rep(16, 3), col=co, legend=c('Vorsø', 'Mågeøerne', 'Stavns Fjord'), bty='n')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
