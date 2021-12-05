# Schaub & Kéry (2022) Integrated Population Models
# Chapter 5 : Introduction to integrated population models
# --------------------------------------------------------

# 5.2 Feeding demographic data into the analysis of a matrix population model
# ===========================================================================

# 5.2.1 Using capture-recapture data in a matrix population model
# ---------------------------------------------------------------

library(IPMbook); library(jagsUI)
data(woodchat5)
str(woodchat5)
# $ ch : num [1:1902, 1:20] 1 1 1 1 1 1 1 1 1 1 ...
# $ age : num [1:1902] 2 2 2 2 2 2 2 2 2 2 ...
# $ repro: num [1:929, 1:3] 6 2 2 5 3 5 3 2 3 2 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:3] "Reproduction" "Year" "Age of mother"
# $ count: num [1:20] 91 119 131 88 139 145 148 116 112 106 ...

marr <- marrayAge(woodchat5$ch, woodchat5$age)

# Bundle data and produce data overview
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), mean.f=c(2.6, 3.6), T=15)
str(jags.data)
# List of 7
# $ marr.j     : num [1:19, 1:20] 8 0 0 0 0 0 0 0 0 0 ...
# $ marr.a     : num [1:19, 1:20] 16 0 0 0 0 0 0 0 0 0 ...
# $ n.occasions: int 20
# $ rel.j      : num [1:19] 51 53 55 65 73 66 61 76 65 75 ...
# $ rel.a      : num [1:19] 36 39 44 61 61 50 43 61 51 53 ...
# $ mean.f     : num [1:2] 2.6 3.6
# $ T          : num 15

# Write JAGS model file
cat(file="model1.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
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

  # Model for population dynamics: a simple matrix population model
  # Define initial stage-specific population sizes
  N[1,1] <- 1
  N[2,1] <- 1
  # Loop over time
  for (t in 1:T){
    # Population projection
    N[1,t+1] <- N[1,t] * mean.f[1] / 2 * mean.sj + N[2,t] * mean.f[2] / 2 * mean.sj
    N[2,t+1] <- (N[1,t] + N[2,t]) * mean.sa
    # Annual growth rate
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])
  }
  lambda <- ann.growth.rate[T]
}
")

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 1), mean.sa=runif(1, 0, 1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "lambda")

# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na)
traceplot(out1)
print(out1, 3)
#             mean    sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# mean.sj    0.304 0.016   0.274   0.304   0.336    FALSE 1 1.001  6000
# mean.sa    0.542 0.014   0.514   0.542   0.569    FALSE 1 1.001  2640
# mean.p     0.604 0.021   0.563   0.604   0.644    FALSE 1 1.001  6000
# lambda     1.018 0.027   0.966   1.018   1.073    FALSE 1 1.001  6000


# 5.2.2 Combining capture-recapture and productivity data in a matrix
#       population model
# -------------------------------------------------------------------

# Bundle data and produce data overview
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
    rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=woodchat5$repro[,1],
    age=woodchat5$repro[,3], T=15)
str(jags.data)
# List of 8
# $ marr.j     : num [1:19, 1:20] 8 0 0 0 0 0 0 0 0 0 ...
# $ marr.a     : num [1:19, 1:20] 16 0 0 0 0 0 0 0 0 0 ...
# $ n.occasions: int 20
# $ rel.j      : num [1:19] 51 53 55 65 73 66 61 76 65 75 ...
# $ rel.a      : num [1:19] 36 39 44 61 61 50 43 61 51 53 ...
# $ J          : num [1:929] 6 2 2 5 3 5 3 2 3 2 ...
# $ age        : num [1:929] 1 1 1 1 1 1 1 1 1 1 ...
# $ T          : num 15

# Write JAGS model file
cat(file="model2.txt", "
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

  # Model for population dynamics: a simple matrix population model
  # Define initial stage-specific population sizes
  N[1,1] <- 1
  N[2,1] <- 1

  # Loop over time
  for (t in 1:T){
    # Population projection
    N[1,t+1] <- N[1,t] * mean.f[1] / 2 * mean.sj + N[2,t] * mean.f[2] / 2 * mean.sj
    N[2,t+1] <- (N[1,t] + N[2,t]) * mean.sa
    # Annual growth rate
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])
  }
  lambda <- ann.growth.rate[T]
}
")

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 1), mean.sa=runif(1, 0, 1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "lambda")

# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out2 <- jags(jags.data, inits, parameters, "model2.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na)
traceplot(out2)
print(out2, 3)
#               mean    sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# mean.sj      0.304 0.016    0.273    0.304    0.337    FALSE 1 1.000  6000
# mean.sa      0.542 0.014    0.514    0.542    0.570    FALSE 1 1.002  1998
# mean.p       0.604 0.021    0.562    0.604    0.643    FALSE 1 1.001  6000
# mean.f[1]    2.677 0.079    2.527    2.675    2.832    FALSE 1 1.001  5055
# mean.f[2]    3.685 0.089    3.517    3.685    3.866    FALSE 1 1.001  2360
# lambda       1.030 0.029    0.973    1.030    1.089    FALSE 1 1.001  5372
