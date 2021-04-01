# Schaub & Kery (2021) Integrated Population Modeling
# Chapter 11 : Woodchat shrike
# ----------------------------
# Code from MS submitted to publisher.

# Run time for test script 4 mins

library(IPMbook) ; library(jagsUI)

# 11.4 Component data likelihoods
# =============================================

library(IPMbook); library(jagsUI)
data(woodchat11)
str(woodchat11)

# List of 6
 # $ ch   : int [1:1079, 1:29] 0 0 0 0 0 0 0 0 0 0 ...
 # $ age  : num [1:1079] 1 1 1 1 1 1 1 1 1 1 ...
 # $ count: num [1:29] 10 13 25 15 17 16 6 16 7 7 ...
 # $ f    : int [1:365] 5 0 6 0 3 6 6 3 0 3 ...
 # $ d    : int [1:365] 0 0 0 0 0 0 0 0 0 1 ...
 # $ year : int [1:365] 1964 1964 1964 1964 1964 1964 1964 1964 1965 1965 ...


# 11.4.1 Population count data (no code)

# 11.4.2 Data on reproduction
# ---------------------------

censored <- woodchat11$f + 1
censored[woodchat11$d==1] <- woodchat11$f[woodchat11$d==1]
woodchat11$f[woodchat11$d==1] <- NA

# 11.4.3 Capture-recapture data
# -----------------------------

marr <- marrayAge(woodchat11$ch, woodchat11$age)
rel.j <- rowSums(marr[,,1])
rel.a <- rowSums(marr[,,2])

recap <- c(1, 2, 3, 4, 5, 6, 7, 8, 19, 9, 10, 11, 12, 19, 19, 19, 19,
    13, 14, 19, 19, 19, 19, 19, 15, 16, 17, 18)

# 11.5 The integrated population model
# ====================================

# Bundle data and produce data overview
jags.data <- list(n.years=dim(marr)[2], marr.j=marr[,,1], marr.a=marr[,,2],
    rel.j=rel.j, rel.a=rel.a, recap=recap, n.recap=max(recap), f=woodchat11$f,
    d=woodchat11$d, censored=censored, year=woodchat11$year-1963,
    C=woodchat11$count, pNinit=dUnif(1, 20))
str(jags.data)

# List of 13
 # $ n.years : int 29
 # $ marr.j  : num [1:28, 1:29] 1 0 0 0 0 0 0 0 0 0 ...
 # $ marr.a  : num [1:28, 1:29] 3 0 0 0 0 0 0 0 0 0 ...
 # $ rel.j   : num [1:28] 27 32 26 32 55 26 22 32 0 25 ...
 # $ rel.a   : num [1:28] 9 14 18 16 15 13 8 17 2 0 ...
 # $ recap   : num [1:28] 1 2 3 4 5 6 7 8 19 9 ...
 # $ n.recap : num 19
 # $ f       : int [1:365] 5 0 6 0 3 6 6 3 0 NA ...
 # $ d       : int [1:365] 0 0 0 0 0 0 0 0 0 1 ...
 # $ censored: num [1:365] 6 1 7 1 4 7 7 4 1 3 ...
 # $ year    : num [1:365] 1 1 1 1 1 1 1 1 2 2 ...
 # $ C       : num [1:29] 10 13 25 15 17 16 6 16 7 7 ...
 # $ pNinit  : num [1:20] 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 ...

# Write JAGS model file
cat(file="model1.txt", "
model {
# Priors and linear models
for (t in 1:(n.years-1)){
  logit.phij[t] ~ dnorm(l.mean.phij, tau.phij)
  phij[t] <- ilogit(logit.phij[t])
  logit.phia[t] ~ dnorm(l.mean.phia, tau.phia)
  phia[t] <- ilogit(logit.phia[t])
  pj[t] <- p.j[recap[t]]
  pa[t] <- p.a[recap[t]]
}
for (t in 1:n.years){
  logit.nu[t] ~ dnorm(l.mean.nu, tau.nu)
  nu[t] <- ilogit(logit.nu[t])
  log.kappa[t] ~ dnorm(l.mean.kappa, tau.kappa)
  kappa[t] <- exp(log.kappa[t])
  log.omega[t] ~ dnorm(l.mean.omega, tau.omega)
  omega[t] <- exp(log.omega[t])
}
for (t in 1:(n.recap-1)){          # Recapture probability
  p.j[t] <- mean.p[1]
  p.a[t] <- mean.p[2]
}
# Fix at 0 recapture probability in years without capture/resighting
p.j[n.recap] <- 0
p.a[n.recap] <- 0

mean.phij ~ dunif(0, 1)
l.mean.phij <- logit(mean.phij)
mean.phia ~ dunif(0, 1)
l.mean.phia <- logit(mean.phia)
mean.nu ~ dunif(0, 1)
l.mean.nu <- logit(mean.nu)
mean.kappa ~ dunif(1, 10)
l.mean.kappa <- log(mean.kappa)
mean.omega ~ dunif(0.01, 30)
l.mean.omega <- log(mean.omega)
for (i in 1:2){
  mean.p[i] ~ dunif(0, 1)
}

sigma.phij ~ dunif(0, 5)
tau.phij <- pow(sigma.phij, -2)
sigma.phia ~ dunif(0, 5)
tau.phia <- pow(sigma.phia, -2)
sigma.nu ~ dunif(0, 5)
tau.nu <- pow(sigma.nu, -2)
sigma.kappa ~ dunif(0, 5)
tau.kappa <- pow(sigma.kappa, -2)
sigma.omega ~ dunif(0, 5)
tau.omega <- pow(sigma.omega, -2)
sigma.f ~ dunif(0, 5)
tau.f <- pow(sigma.f, -2)

# Population count data (state-space model)
# Model for the initial population size: uniform priors
R[1] ~ dcat(pNinit)      # Local recruits
S[1] ~ dcat(pNinit)      # Surviving adults
I[1] ~ dpois(omega[1])   # Immigrants

# Process model over time: our model of population dynamics
for (t in 2:n.years){
  R[t] ~ dpois(nu[t-1] * kappa[t-1]/2 * phij[t-1] * N[t-1]) # Local recruits
  S[t] ~ dbin(phia[t-1], N[t-1])     # Surviving adults
  I[t] ~ dpois(omega[t])             # Immigrants
}

# Observation model
for (t in 1:n.years){
  N[t] <- S[t] + R[t] + I[t]         # Total number of breeding females
  C[t] ~ dpois(N[t])

  # GoF for population count data: mean absolute percentage error
  C.pred[t] ~ dpois(N[t])
  disc.C[t] <- pow(((C[t] - N[t]) / C[t]) * ((C[t] - N[t]) / (C[t] + 0.001)), 0.5)   # Add a small number to avoid potential division by 0
  discN.C[t] <- pow(((C.pred[t] - N[t]) / (C.pred[t] + 0.001)) * ((C.pred[t] - N[t]) / (C.pred[t] + 0.001)), 0.5)
}
fit.C <- 100 / n.years * sum(disc.C)
fitN.C <- 100 / n.years * sum(discN.C)

# Productivity data (zero-inflated normal with right censoring)
for (i in 1:length(f)){
  z[i] ~ dbern(nu[year[i]])
  f[i] ~ dnorm(z[i] * kappa[year[i]], tau.f)
  delta[i] <- step(f[i] - censored[i])
  d[i] ~ dbern(z[i] * delta[i])

  # GoF for productivity data: deviance
  z.pred[i] ~ dbern(nu[year[i]])
  f.pred[i] ~ dnorm(z.pred[i] * kappa[year[i]], tau.f)
  f.exp[i] <- kappa[year[i]] * nu[year[i]]

  dev[i] <- (f[i]  - f.exp[i]) * (f[i]  - f.exp[i])
  devN[i] <- (f.pred[i] - f.exp[i]) * (f.pred[i] - f.exp[i])
}
fit.f <- 2 * sum(dev)
fitN.f <- 2 * sum(devN)

# Capture-recapture data (CJS model with multinomial likelihood)
# Define the multinomial likelihood
for (t in 1:(n.years-1)){
  marr.j[t,1:n.years] ~ dmulti(pr.j[t,], rel.j[t])
  marr.a[t,1:n.years] ~ dmulti(pr.a[t,], rel.a[t])
}
# Define the cell probabilities of the m-arrays
# Main diagonal
for (t in 1:(n.years-1)){
  qj[t] <- 1-pj[t]
  qa[t] <- 1-pa[t]
  pr.j[t,t] <- phij[t] * pj[t]
  pr.a[t,t] <- phia[t] * pa[t]
  # Above main diagonal
  for (j in (t+1):(n.years-1)){
    pr.j[t,j] <- phij[t] * prod(phia[(t+1):j]) * qj[t] * prod(qa[t:(j-1)]) * pa[j] / qa[t]
    pr.a[t,j] <- prod(phia[t:j]) * prod(qa[t:(j-1)]) * pa[j]
  } #j
  # Below main diagonal
  for (j in 1:(t-1)){
    pr.j[t,j] <- 0
    pr.a[t,j] <- 0
  } #j
} #t
# Last column: probability of non-recapture
for (t in 1:(n.years-1)){
  pr.j[t,n.years] <- 1-sum(pr.j[t,1:(n.years-1)])
  pr.a[t,n.years] <- 1-sum(pr.a[t,1:(n.years-1)])
}

# GoF for cap.-recap.data: Freeman-Tukey test statistics
for (t in 1:(n.years-1)){
  # Simulated m-arrays
  marr.j.pred[t,1:n.years] ~ dmulti(pr.j[t,], rel.j[t])
  marr.a.pred[t,1:n.years] ~ dmulti(pr.a[t,], rel.a[t])

  # Expected values and test statistics
  for (j in 1:n.years){
    marr.jE[t,j] <- pr.j[t,j] * rel.j[t]
    E.org.j[t,j] <- pow((pow(marr.j[t,j], 0.5) - pow(marr.jE[t,j], 0.5)), 2)
    E.new.j[t,j] <- pow((pow(marr.j.pred[t,j], 0.5) - pow(marr.jE[t,j], 0.5)), 2)
    marr.aE[t,j] <- pr.a[t,j] * rel.a[t]
    E.org.a[t,j] <- pow((pow(marr.a[t,j], 0.5) - pow(marr.aE[t,j], 0.5)), 2)
    E.new.a[t,j] <- pow((pow(marr.a.pred[t,j], 0.5) - pow(marr.aE[t,j], 0.5)), 2)
  } #j
} #t

fit.CJS <- sum(E.org.j) + sum(E.org.a)
fitN.CJS <- sum(E.new.j) + sum(E.new.a)

# Compute derived parameters: total productivity and immigration rate
for (t in 1:n.years){
  rho[t] <- nu[t] * kappa[t]
}
mean.rho <- mean.nu * mean.kappa
for (t in 1:(n.years-1)){
  omega.rate[t] <- I[t+1] / N[t]
}
}
")

# Initial values
init.Juv <- censored+1
init.Juv[woodchat11$d==0] <- NA
s <- rep(1, length(woodchat11$f))
s[woodchat11$f==0] <- 0
inits <- function() {list(z=s, f=init.Juv)}

# Parameters monitored
parameters <- c("mean.phij", "mean.phia", "mean.nu", "mean.kappa", "mean.rho",
    "mean.omega", "sigma.phij", "sigma.phia", "sigma.nu", "sigma.kappa",
    "sigma.omega", "sigma.f", "mean.p", "fit.C", "fitN.C", "fit.f", "fitN.f",
    "fit.CJS", "fitN.CJS", "phij", "phia", "nu", "kappa", "omega", "rho",
    "omega.rate", "R", "S", "I", "N")

# MCMC settings
ni <- 30000; nb <- 10000; nc <- 3; nt <- 10; na <- 1000

# Call JAGS from R (ART 3 min) and check convergence
out1 <- jags(jags.data, inits, parameters, "model1.txt",
    n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)
par(mfrow=c(3, 3)); traceplot(out1)


# 11.6 Results
# ============

print(out1, 3)

                   # mean       sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# mean.phij         0.052    0.024    0.016    0.048    0.107    FALSE 1 1.017   178
# mean.phia         0.411    0.068    0.287    0.407    0.556    FALSE 1 1.005  2004
# mean.nu           0.742    0.037    0.666    0.742    0.812    FALSE 1 1.000  6000
# mean.kappa        4.757    0.083    4.593    4.757    4.923    FALSE 1 1.000  6000
# mean.rho          3.528    0.184    3.151    3.534    3.875    FALSE 1 1.000  6000
# mean.omega        6.817    1.424    4.370    6.683    9.600    FALSE 1 1.015  2040
# sigma.phij        0.566    0.461    0.025    0.462    1.699    FALSE 1 1.010   287
# sigma.phia        0.384    0.284    0.020    0.330    1.054    FALSE 1 1.002  1012
# sigma.nu          0.705    0.220    0.293    0.696    1.156    FALSE 1 1.003   705
# sigma.kappa       0.056    0.021    0.015    0.056    0.097    FALSE 1 1.001  3637
# sigma.omega       0.276    0.198    0.019    0.242    0.731    FALSE 1 1.004   614
# sigma.f           0.879    0.037    0.809    0.877    0.956    FALSE 1 1.000  4028
# mean.p[1]         0.147    0.082    0.037    0.130    0.355    FALSE 1 1.002  1047
# mean.p[2]         0.462    0.104    0.277    0.457    0.681    FALSE 1 1.000  6000
# fit.C            22.401    3.506   16.117   22.290   29.613    FALSE 1 1.000  6000
# fitN.C           74.001 1053.304   15.973   24.282   39.925    FALSE 1 1.005  6000
# fit.f          3378.453   92.348 3221.118 3370.536 3586.062    FALSE 1 1.002  1319
# fitN.f         3593.314  311.470 2992.793 3589.004 4213.025    FALSE 1 1.000  5265
# fit.CJS          23.864    3.025   18.618   23.685   30.302    FALSE 1 1.002  1333
# fitN.CJS         22.754    3.927   15.495   22.650   30.811    FALSE 1 1.001  1443
# phij[1]           0.066    0.044    0.015    0.056    0.180    FALSE 1 1.007   619
# [... output truncated ...]
# phij[28]          0.056    0.041    0.006    0.048    0.159    FALSE 1 1.009   484
# phia[1]           0.483    0.127    0.291    0.460    0.781    FALSE 1 1.001  5620
# [... output truncated ...]
# phia[28]          0.404    0.109    0.197    0.398    0.650    FALSE 1 1.001  1235
# nu[1]             0.740    0.099    0.524    0.749    0.903    FALSE 1 1.000  6000
# [... output truncated ...]
# nu[29]            0.779    0.091    0.578    0.788    0.928    FALSE 1 1.000  3355
# kappa[1]          4.792    0.221    4.363    4.786    5.251    FALSE 1 1.000  6000
# [... output truncated ...]
# kappa[29]         4.602    0.222    4.146    4.615    5.008    FALSE 1 1.000  6000
# omega[1]          6.433    2.196    2.282    6.368   10.739    FALSE 1 1.005  2730
# [... output truncated ...]
# omega[29]         6.667    2.185    2.576    6.591   11.282    FALSE 1 1.003  6000
# rho[1]            3.546    0.500    2.480    3.586    4.409    FALSE 1 1.000  6000
# [... output truncated ...]
# rho[29]           3.583    0.438    2.640    3.611    4.338    FALSE 1 1.001  2784
# omega.rate[1]     0.567    0.314    0.115    0.500    1.333    FALSE 1 1.001  2628
# [... output truncated ...]
# omega.rate[28]    0.588    0.312    0.125    0.545    1.333    FALSE 1 1.000  6000
# R[1]              4.192    2.921    1.000    4.000   11.000    FALSE 1 1.002  1537
# [... output truncated ...]
# R[29]             0.843    1.087    0.000    1.000    4.000     TRUE 1 1.002  1436
# S[1]              4.149    2.912    1.000    3.000   11.000    FALSE 1 1.000  6000
# [... output truncated ...]
# S[29]             4.214    1.996    1.000    4.000    9.000    FALSE 1 1.000  5499
# I[1]              5.138    2.590    1.000    5.000   11.000    FALSE 1 1.000  6000
# [... output truncated ...]
# I[29]             6.007    2.632    2.000    6.000   12.000    FALSE 1 1.001  6000
# N[1]             13.480    3.209    8.000   13.000   20.000    FALSE 1 1.001  4489
# [... output truncated ...]
# N[29]            11.063    2.641    6.000   11.000   17.000    FALSE 1 1.000  6000

# Compute CRI of juvenile survival and emigration
sj <- 2 * (1 - out1$sims.list$mean.phia) / out1$sims.list$mean.nu / out1$sims.list$mean.kappa
mean(sj); quantile(sj, c(0.025, 0.975))

xsi <- 1 - out1$sims.list$mean.phij / sj
mean(xsi); quantile(xsi, c(0.025, 0.975))

# ~~~ save output ~~~
save(out1, file="IPM_11_output.RData")
# ~~~~~~~~~~~~~~~~~~~
