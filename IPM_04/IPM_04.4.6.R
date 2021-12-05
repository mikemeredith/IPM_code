# Schaub & Kéry (2022) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------

# 4.4 Models for productivity surveys
# ===================================

# 4.4.6 Nest survival models
# --------------------------

library(IPMbook); library(jagsUI)
data(wryneck)
str(wryneck)
# 'data.frame': 181 obs. of 5 variables:
# $ f : int 49 15 23 23 47 11 23 15 65 2 ...
# $ j : int 53 32 37 40 66 30 37 30 71 18 ...
# $ k : int 53 32 37 40 66 30 37 30 72 18 ...
# $ x : int 1 1 1 1 1 1 1 1 0 1 ...
# $ age: int 16 2 4 2 2 2 5 2 5 2 ...

tail(wryneck)
#      f  j  k x age
# 176  5 25 25 1   1
# 177 53 56 60 0   1
# 178  9 29 29 1   1
# 179 26 46 46 1   1
# 180 55 75 75 1   1
# 181 20 35 35 1   6

# Identify failed broods
fail <- which(wryneck$x==0)

# Create encounter histories
y <- matrix(NA, nrow=length(wryneck$f), ncol=max(wryneck$k))
for (i in 1:length(wryneck$f)){
  y[i,wryneck$f[i]] <- 1
  y[i,wryneck$j[i]] <- 1
}
for (i in 1:length(fail)){
  y[fail[i],wryneck$k[fail[i]]] <- 0
}

y[176,]
# [1] ... NA 1 NA NA ... NA NA 1 NA NA ...

y[177,]
# [37] ... NA NA 1 NA NA 1 NA NA NA 0 NA NA ...

# Bundle data
jags.data <- with(wryneck, list(y=y, f=f, k=k, n.nest=nrow(y), T=21, age=age))
str(jags.data)
# List of 6
# $ y     : num [1:181, 1:81] NA NA NA NA NA NA NA NA NA NA ...
# $ f     : int [1:181] 49 15 23 23 47 11 23 15 65 2 ...
# $ k     : int [1:181] 53 32 37 40 66 30 37 30 72 18 ...
# $ n.nest: int 181
# $ T     : num 21
# $ age   : int [1:181] 16 2 4 2 2 2 5 2 5 2 ...

# Write JAGS model file
cat(file="model13.txt", "
model {
  # Priors and linear models
  for (i in 1:n.nest){
    for (t in f[i]:(k[i]-1)){
      phi[i,t] <- phia[age[i] + t - f[i]]
    } #t
  } #i

  for (a in 1:T){
    phia[a] <- ilogit(alpha + beta * a)
  }
  alpha ~ dnorm(0, 0.001)
  beta ~ dnorm(0, 0.001)

  # Likelihood
  for (i in 1:n.nest){
    for (t in (f[i]+1):k[i]){
      y[i,t] ~ dbern(phi[i,t-1] * y[i,t-1])
    } #t
  } #i

  # Derived parameter: nest success
  nu <- prod(phia[1:T])
}
")

# Initial values
inits <- function(){list(alpha=runif(1, 4, 5), beta=runif(1, 0, 0.1))}

# Parameters monitored
parameters <- c("phia", "nu", "alpha", "beta")

# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000

# Call JAGS from R (ART 1 min) and check convergence
out16 <- jags(jags.data, inits, parameters, "model13.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out16) # Not shown
print(out16, 3)
#             mean     sd   2.5%     50%   97.5% overlap0     f  Rhat n.eff
# phia[1]    0.971  0.008  0.953   0.972   0.985    FALSE 1.000 1.001  1650
# phia[2]    0.973  0.007  0.958   0.974   0.985    FALSE 1.000 1.001  1633
# [... output truncated ...]
# phia[20]   0.993  0.003  0.987   0.993   0.997    FALSE 1.000 1.001  5757
# phia[21]   0.993  0.003  0.987   0.994   0.997    FALSE 1.000 1.001  4989
# nu         0.733  0.036  0.659   0.733   0.800    FALSE 1.000 1.000  3483
# alpha      3.479  0.319  2.871   3.470   4.137    FALSE 1.000 1.002  1361
# beta       0.078  0.032  0.016   0.078   0.141    FALSE 0.995 1.001  2152

# ~~~ code for figure 4.13 ~~~~
plot(out16$mean$phia, type="b", pch=16, ylim=range(c(out16$q2.5$phia, 1)),
    ylab="Daily nest survival probability", xlab="Nestling age", axes=FALSE)
segments(1:21, out16$q2.5$phia, 1:21, out16$q97.5$phia)
axis(1, at=1:21, labels=NA, tcl=-0.25)
axis(1, at=c(5,10,15,20), labels=c("5","10","15","20"), tcl=-0.5)
axis(2, las=1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
