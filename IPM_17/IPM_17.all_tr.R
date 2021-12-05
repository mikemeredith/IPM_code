# Schaub & Kéry (2022) Integrated Population Models
# Chapter 17 : Elk
# ----------------

# Run time for test script 26 mins, full run 65 mins

# 17.4 Component data likelihoods
# =============================================

library(IPMbook); library(jagsUI)
data(elk)
str(elk)
# List of 3
# $ C: num [1:17, 1:6] 21 36 15 8 7 8 6 5 3 3 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:17] "1" "2" "3" "4" ...
# .. ..$ : chr [1:6] "1988" "1989" "1990" "1991" ...
# $ H: num [1:2, 1:6] 275 143 290 154 211 211 360 272 201 201 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:2] "Hunters surveyed" "Registered harvest"
# .. ..$ : chr [1:6] "1988" "1989" "1990" "1991" ...
# $ R: num [1:3, 1:6] 1 2 10 4 1 20 0 1 24 0 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:3] "hunted" "died naturally" "survived"
# .. ..$ : chr [1:6] "1988" "1989" "1990" "1991" ...


# 17.4.1 Age-at-harvest data (no code)
# 17.4.2 Hunter survey data (no code)
# 17.4.3 Radio tracking data (no code)


# 17.5 The integrated population model
# ====================================

s <- 0.8                                          # Guess of annual survival
f <- 0.25                                         # Guess of recruitment

# Create matrix population model (transition matrix A)
A <- matrix(0, ncol=17, nrow=17)
A[1,2:17] <- f
for (a in 2:17){
  A[a,a-1] <- s
}

# Compute stable age distribution (right eigenvector of A)
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])

# Population size in first age class of first year
r <- 0.5                                          # Guess of reporting rate
h <- 0.1                                          # Guess of hunting mortality
N1 <- elk$C[1,1] / (h * r)

# Compute age-specific population sizes in first year
n <- N1 * revec / revec[1]

# Bundle data and produce data overview
jags.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R, total=colSums(R), n.years=ncol(C),
    n.age=nrow(C), n=n))
str(jags.data)
# List of 8
# $ C      : num [1:17, 1:6] 21 36 15 8 7 8 6 5 3 3 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:17] "1" "2" "3" "4" ...
# .. ..$ : chr [1:6] "1988" "1989" "1990" "1991" ...
# $ a      : Named num [1:6] 275 290 211 360 201 325
# ..- attr(*, "names")= chr [1:6] "1988" "1989" "1990" "1991" ...
# $ b      : Named num [1:6] 143 154 211 272 201 155
# ..- attr(*, "names")= chr [1:6] "1988" "1989" "1990" "1991" ...
# $ R      : num [1:3, 1:6] 1 2 10 4 1 20 0 1 24 0 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:3] "hunted" "died naturally" "survived"
# .. ..$ : chr [1:6] "1988" "1989" "1990" "1991" ...
# $ total  : Named num [1:6] 13 25 25 16 25 23
# ..- attr(*, "names")= chr [1:6] "1988" "1989" "1990" "1991" ...
# $ n.years: int 6
# $ n.age  : int 17
# $ n      : num [1:17] 420 338 272 218 176 ...

# Write JAGS model file
cat(file = "model1.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.years-1)){
    f[t] ~ dunif(0, 1)
  }
  for (t in 1:n.years){
    s[t] ~ dunif(0, 1)
    h[t] ~ dunif(0, 1)
    r[t] ~ dunif(0, 1)
  }

  # Age-at-harvest data (state-space model)
  # Model for the initial population size: Poisson priors
  for (a in 1:n.age){
    N[a,1] ~ dpois(n[a])
  }

  # Process model over time: our model of population dynamics
  for (t in 1:(n.years-1)){
    N[1,t+1] ~ dpois((Ntot[t] - N[1,t]) * f[t])
    for (a in 2:n.age){
      N[a,t+1] ~ dbin((1-h[t]) * s[t], N[a-1,t])
    } #a
  } #t

  # Derived quantity: total year-specific population size
  for (t in 1:n.years){
    Ntot[t] <- sum(N[,t])
  }

  # Observation model
  for (t in 1:n.years){
    for (a in 1:n.age){
      C[a,t] ~ dbin(h[t] * r[t], N[a,t])
    } #a
  } #t

  # Hunter survey data (logistic regression model)
  for (t in 1:n.years){
    b[t] ~ dbin(r[t], a[t])
  }

  # Radio tracking data (multinomial model)
  for (t in 1:n.years){
    R[,t] ~ dmulti(prt[,t], total[t])
    prt[1,t] <- h[t]
    prt[2,t] <- (1-h[t]) * (1-s[t])
    prt[3,t] <- (1-h[t]) * s[t]
  }
}
")

# Initial values
inits <- function() {list(s=runif(jags.data$n.years, 0.8, 1), h=runif(jags.data$n.years,
    0.05, 0.15))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
ni <- 250000; nb <- 50000; nc <- 3; nt <- 20; na <- 5000

# Call JAGS from R (ART 4 min) and check convergence
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1)


# 17.6 Results on elk population dynamics
# =======================================

print(out1, 3)

#              mean      sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# h[1]        0.116   0.012    0.094    0.115    0.141    FALSE 1 1.000 30000
# h[2]        0.145   0.018    0.114    0.144    0.184    FALSE 1 1.001  6631
# h[3]        0.086   0.011    0.067    0.085    0.111    FALSE 1 1.001  3978
# h[4]        0.085   0.013    0.063    0.084    0.115    FALSE 1 1.002  1590
# h[5]        0.080   0.013    0.058    0.079    0.110    FALSE 1 1.001  1802
# h[6]        0.182   0.036    0.124    0.178    0.265    FALSE 1 1.003  1062
# s[1]        0.851   0.068    0.701    0.857    0.963    FALSE 1 1.001  2988
# s[2]        0.929   0.046    0.817    0.938    0.991    FALSE 1 1.000 10767
# s[3]        0.934   0.044    0.824    0.943    0.991    FALSE 1 1.001  7929
# s[4]        0.951   0.046    0.828    0.965    0.999    FALSE 1 1.001  4161
# s[5]        0.817   0.071    0.662    0.823    0.935    FALSE 1 1.001  1967
# s[6]        0.792   0.082    0.611    0.800    0.926    FALSE 1 1.000 30000
# f[1]        0.259   0.035    0.193    0.258    0.331    FALSE 1 1.000  4754
# f[2]        0.300   0.038    0.231    0.298    0.378    FALSE 1 1.000 30000
# f[3]        0.316   0.043    0.240    0.314    0.407    FALSE 1 1.000  4837
# f[4]        0.207   0.034    0.145    0.205    0.280    FALSE 1 1.000 29810
# f[5]        0.140   0.034    0.083    0.137    0.215    FALSE 1 1.001  2521
# r[1]        0.520   0.030    0.462    0.520    0.579    FALSE 1 1.000 30000
# r[2]        0.530   0.029    0.473    0.530    0.587    FALSE 1 1.000 14354
# r[3]        0.995   0.005    0.983    0.997    1.000    FALSE 1 1.000 15095
# r[4]        0.755   0.023    0.709    0.755    0.797    FALSE 1 1.000 30000
# r[5]        0.995   0.005    0.982    0.997    1.000    FALSE 1 1.001 10259
# r[6]        0.482   0.027    0.429    0.482    0.536    FALSE 1 1.000 22795
# Ntot[1]  2096.276  45.462 2007.000 2096.000 2187.000    FALSE 1 1.000 23146
# Ntot[2]  1997.991 170.065 1641.975 2007.000 2302.025    FALSE 1 1.001  2137
# Ntot[3]  2048.117 222.540 1608.000 2051.000 2478.000    FALSE 1 1.001  1819
# Ntot[4]  2239.809 285.800 1693.000 2236.000 2815.000    FALSE 1 1.002  1255
# Ntot[5]  2299.903 336.905 1675.000 2291.000 2995.025    FALSE 1 1.001  1903
# Ntot[6]  1991.210 355.567 1352.000 1973.000 2742.000    FALSE 1 1.002   964
# N[1,1]    428.762  18.718  393.000  429.000  466.000    FALSE 1 1.000  9726
# N[2,1]    355.155  16.736  323.000  355.000  388.000    FALSE 1 1.000  8812
# [...output truncated...]
# N[16,6]    15.846   4.218    8.000   16.000   25.000    FALSE 1 1.001  2284
# N[17,6]    12.726   3.854    6.000   12.000   21.000    FALSE 1 1.001  3377


# Calculation of annual survival (s*)
s.star <- (1-out1$sims.list$h) * out1$sims.list$s

# ~~~~ Code for Fig. 17.3 ~~~~
library(scales)
cl <- viridis_pal(option='E')(20)[c(18,2)]
qu <- function(x) quantile(x, c(0.025, 0.975))
op <- par(mfrow=c(2,1), mar=c(3.5, 5, 1, 1))
d <- 0.1
plot(y=out1$mean$h, x=(1:6)+d, ylim=c(0,0.4), xlim=c(1-d, 6+d), type="b",
    pch=16, axes=FALSE, ylab="Probability", xlab=NA, col=cl[2])
segments((1:6)+d, out1$q2.5$h, (1:6)+d, out1$q97.5$h, col=cl[2])
axis(1, at=1:6, labels = 1988:1993)
axis(2, las=1)

points(y=1-out1$mean$s, x=(1:6)-d, type="b", pch=16, col=cl[1])
segments((1:6)-d, 1-out1$q2.5$s, (1:6)-d, 1-out1$q97.5$s, col=cl[1])
legend('topleft', pch=rep(16,2), col=rev(cl),
    legend=c('Hunting mortality', 'Background mortality'), bty='n')

plot(out1$mean$f, x=2:6, ylim=c(0,1.2), xlim=c(1-d, 6+d), type="b",
    pch=16, axes=FALSE, ylab="Probability", xlab=NA, col=cl[1])
segments(2:6, out1$q2.5$f, 2:6, out1$q97.5$f, col=cl[1])
axis(1, at=1:6, labels = 1988:1993)
axis(2, las=1, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1))
points(y=apply(s.star, 2, mean), x=1:6, type="b", pch=16, col=cl[2])
segments(1:6, apply(s.star, 2, qu)[1,], 1:6, apply(s.star, 2, qu)[2,], col=cl[2])
legend('topleft', pch=c(16,16), col=rev(cl),
    legend=c('Annual survival', 'Recruitment'), bty='n')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 17.4 ~~~~
library(scales)
cl <- viridis_pal(option='E')(20)[c(18,11)]
qu <- function(x) quantile(x, c(0.025, 0.975))

# Calculate the number of harvested elk
H <- out1$sims.list$Ntot * out1$sims.list$h

op <- par(mar=c(3.5, 5, 1, 1), las=1)
z <- cbind(out1$sims.list$Ntot[,1], H[,1], out1$sims.list$Ntot[,2], H[,2],  # ~~~ FIXME
    out1$sims.list$Ntot[,3], H[,3], out1$sims.list$Ntot[,4], H[,4],
    out1$sims.list$Ntot[,5], H[,5], out1$sims.list$Ntot[,6], H[,6])
a <- barplot(apply(z, 2, mean), space = rep(c(0.5,0.1), 6), ylim = c(0, 3000),
    ylab="Numbers", col=rep(cl,6), border=NA)
axis(1, at = c(mean(a[1:2]), mean(a[3:4]), mean(a[5:6]), mean(a[7:8]),
    mean(a[9:10]), mean(a[11:12])), labels = 1988:1993, lwd=0)
segments(a, apply(z, 2, qu)[1,], a, apply(z, 2, qu)[2,], col='black')
legend('topleft', pch=c(15,15), legend=c('Pre-hunting population', 'Harvested'),
    bty='n', col=cl, pt.cex=1.75)
segments(a[1]-0.55, 0, a[12]+0.55, 0)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 17.7 Prior sensitivity analysis
# ===============================

# ~~~~ code for analysis with different priors ~~~~
# Prior Poisson 2 (P2): h = 0.075
s <- 0.8   # guess of annual survival
f <- 0.25  # guess of recruitment

# Create matrix population model
A <- matrix(0, ncol=17, nrow=17)
A[1,2:17] <- f
for (a in 2:17){
  A[a,a-1] <- s
}

# Compute right eigenvector
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])

# Population size in first age class of first year
r <- 0.5     # guess of reporting rate
h <- 0.075   # guess of hunting mortality
N1 <- elk$C[1,1] / (h * r)

# Compute age-dependent population size in first year
n <- N1 * revec / revec[1]

# Bundle data and produce data overview
jags.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R, total=colSums(R),
    n.years=ncol(C), n.age=nrow(C), n=n))

# Initial values
inits <- function() {list(s=runif(jags.data$n.years, 0.8, 1),
    h=runif(jags.data$n.years, 0.05, 0.15))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
ni <- 250000; nt <- 20; nb <- 50000; nc <- 3; na <- 5000

# Call JAGS from R (ART 4 min) and check convergence
out2 <- jags(jags.data, inits, parameters, "model1.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = na, parallel = TRUE)


# Prior Poisson 3 (P3): h = 0.125
s <- 0.8   # guess of annual survival
f <- 0.25  # guess of recruitment

# Create matrix population model
A <- matrix(0, ncol=17, nrow=17)
A[1,2:17] <- f
for (a in 2:17){
  A[a,a-1] <- s
}

# Compute right eigenvector
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])

# Population size in first age class of first year
r <- 0.5     # guess of reporting rate
h <- 0.125   # guess of hunting mortality
N1 <- elk$C[1,1] / (h * r)

# Compute age-dependent population size in first year
n <- N1 * revec / revec[1]

# Bundle data and produce data overview
jags.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R, total=colSums(R),
    n.years=ncol(C), n.age=nrow(C), n=n))

# Initial values
inits <- function() {list(s=runif(jags.data$n.years, 0.8, 1), h=runif(jags.data$n.years, 0.05, 0.15))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
ni <- 250000; nt <- 20; nb <- 50000; nc <- 3; na <- 5000

# Call JAGS from R (ART 4 min) and check convergence
out3 <- jags(jags.data, inits, parameters, "model1.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = na, parallel = TRUE)


# Discrete uniform prior (U)
s <- 0.8   # guess of annual survival
f <- 0.25  # guess of recruitment

# Create matrix population model
A <- matrix(0, ncol=17, nrow=17)
A[1,2:17] <- f
for (a in 2:17){
  A[a,a-1] <- s
}

# Compute right eigenvector
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])

# Population size in first age class of first year
r <- 0.5     # guess of reporting rate
h <- 0.1     # guess of hunting mortality
N1 <- elk$C[1,1] / (h * r)

# Compute age-dependent population size in first year
n <- N1 * revec / revec[1]
lo <- 0.5 * n
up <- 1.5 * n

# Create a matrix with the values for the categorical distribution
pinit <- dUnif(lower=lo, upper=up)

# Bundle data
jags.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R, total=colSums(R),
    n.years=ncol(C), n.age=nrow(C), p=pinit, up=round(up)))


# Write JAGS model file
cat(file = "model2.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.years-1)){
    f[t] ~ dunif(0, 1)
  }
  for (t in 1:n.years){
    s[t] ~ dunif(0, 1)
    h[t] ~ dunif(0, 1)
    r[t] ~ dunif(0, 1)
  }

  # Age-at-harvest data (state-space model)
  # Model for the initial population size: discrete uniform priors
  for (i in 1:17) {
    N[i,1] ~ dcat(p[i, 1:up[i]])
  }
  # Process model over time: our model of population dynamics
  for (t in 1:(n.years-1)){
    N[1,t+1] ~ dpois((Ntot[t] - N[1,t]) * f[t])
    for (a in 2:n.age){
      N[a,t+1] ~ dbin((1-h[t]) * s[t], N[a-1,t])
    } #a
  } #t

  # Derived quantity: total year-specific population size
  for (t in 1:n.years){
    Ntot[t] <- sum(N[,t])
  }

  # Observation model
  for (t in 1:n.years){
    for (a in 1:n.age){
      C[a,t] ~ dbin(h[t] * r[t], N[a,t])
    } #a
  } #t

  # Hunter survey data (logistic regression model)
  for (t in 1:n.years){
    b[t] ~ dbin(r[t], a[t])
  }

  # Radio tracking data (multinomial model)
  for (t in 1:n.years){
    R[,t] ~ dmulti(prt[,t], total[t])
    prt[1,t] <- h[t]
    prt[2,t] <- (1-h[t]) * (1-s[t])
    prt[3,t] <- (1-h[t]) * s[t]
  }
}
")

# Initial values
inits <- function() {list(s=runif(jags.data$n.years, 0.8, 1),
    h=runif(jags.data$n.years, 0.05, 0.15))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
# ni <- 250000; nt <- 20; nb <- 50000; nc <- 3; na <- 5000
ni <- 25000; nt <- 2; nb <- 5000; nc <- 3; na <- 500  # ~~~ for testing

# Call JAGS from R (ART 16 min) and check convergence
out4 <- jags(jags.data, inits, parameters, "model2.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = na, parallel = TRUE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 17.5 ~~~~
library(scales)
cl <- viridis_pal(option='E')(20)[c(18,11,5,1)]
op <- par(mfrow=c(2,1), las=1, mar=c(3,4,1,1))
d <- 0.1
plot(x=(1:6)-1.5*d, y=out1$mean$Ntot, ylim = c(0, 4000), xlim=c(1-1.5*d, 6+1.5*d), pch=16,
    axes=FALSE, xlab=NA, ylab=expression(paste("Population size (", italic(N)[tot],")")), col=cl[1])
segments((1:6)-1.5*d, out1$q2.5$Ntot, (1:6)-1.5*d, out1$q97.5$Ntot, col=cl[1])
points(x=(1:6)-0.5*d, y=out2$mean$Ntot, pch=16, col=cl[2])
segments((1:6)-0.5*d, out2$q2.5$Ntot, (1:6)-0.5*d, out2$q97.5$Ntot, col=cl[2])
points(x=(1:6)+0.5*d, y=out3$mean$Ntot, pch=16, col=cl[3])
segments((1:6)+0.5*d, out3$q2.5$Ntot, (1:6)+0.5*d, out3$q97.5$Ntot, col=cl[3])
points(x=(1:6)+1.5*d, y=out4$mean$Ntot, pch=16, col=cl[4])
segments((1:6)+1.5*d, out4$q2.5$Ntot, (1:6)+1.5*d, out4$q97.5$Ntot, col=cl[4])
axis(1, at=1:6, labels=1988:1993)
axis(2)

plot(x=(1:6)-1.5*d, y=out1$mean$h, ylim = c(0, 0.31), xlim=c(1-1.5*d, 6+1.5*d), pch=16,
    axes=FALSE, xlab=NA, ylab=expression(paste("Hunting mortality (", italic(h), ")")), col=cl[1])
segments((1:6)-1.5*d, out1$q2.5$h, (1:6)-1.5*d, out1$q97.5$h, col=cl[1])
points(x=(1:6)-0.5*d, y=out2$mean$h, pch=16, col=cl[2])
segments((1:6)-0.5*d, out2$q2.5$h, (1:6)-0.5*d, out2$q97.5$h, col=cl[2])
points(x=(1:6)+0.5*d, y=out3$mean$h, pch=16, col=cl[3])
segments((1:6)+0.5*d, out3$q2.5$h, (1:6)+0.5*d, out3$q97.5$h, col=cl[3])
points(x=(1:6)+1.5*d, y=out4$mean$h, pch=16, col=cl[4])
segments((1:6)+1.5*d, out4$q2.5$h, (1:6)+1.5*d, out4$q97.5$h, col=cl[4])
axis(1, at=1:6, labels=1988:1993)
axis(2)
legend(x=0.75, y=0.33, pch=rep(16,2), col=cl[1:2],
    legend=c(expression('P'[1]),expression('P'[2])), bty='n')
legend(x=1.5, y=0.33, , inset = 0.1, pch=rep(16,2), col=cl[3:4],
    legend=c(expression('P'[3]), 'U'), bty='n')
par(op)

# Fig. 17.6
library(scales)
cl <- viridis_pal(option='E')(20)[c(18,11,5,1)]
qu <- function(x) quantile(x, c(0.025, 0.975))

# Calculate annual population growth rates
lam1 <- lam2 <- lam3 <- lam4 <-
    matrix(NA, ncol=ncol(out1$sims.list$Ntot)-1, nrow=nrow(out1$sims.list$Ntot))
for (t in 1:5){
  lam1[,t] <- out1$sims.list$Ntot[,t+1] / out1$sims.list$Ntot[,t]
  lam2[,t] <- out2$sims.list$Ntot[,t+1] / out2$sims.list$Ntot[,t]
  lam3[,t] <- out3$sims.list$Ntot[,t+1] / out3$sims.list$Ntot[,t]
  lam4[,t] <- out4$sims.list$Ntot[,t+1] / out4$sims.list$Ntot[,t]
}

op <- par(las = 1, mar = c(3,4,1,1))
d <- 0.1
xa <- c(1.5, 2.5, 3.5, 4.5, 5.5)
plot(x=xa-1.5*d, y=apply(lam1, 2, mean), type='b', pch=16, xlim=c(1,6), ylim=c(0.6, 1.2),
    col=cl[1], axes=FALSE, xlab=NA, ylab='Population growth rate')
segments(xa-1.5*d, apply(lam1, 2, qu)[1,], xa-1.5*d, apply(lam1, 2, qu)[2,], col=cl[1])
points(x=xa-0.5*d, y=apply(lam2, 2, mean), type='b', pch=16, col=cl[2])
segments(xa-0.5*d, apply(lam2, 2, qu)[1,], xa-0.5*d, apply(lam2, 2, qu)[2,], col=cl[2])
points(x=xa+0.5*d, y=apply(lam3, 2, mean), type='b', pch=16, col=cl[3])
segments(xa+0.5*d, apply(lam3, 2, qu)[1,], xa+0.5*d, apply(lam3, 2, qu)[2,], col=cl[3])
points(x=xa+1.5*d, y=apply(lam4, 2, mean), type='b', pch=16, col=cl[4])
segments(xa+1.5*d, apply(lam4, 2, qu)[1,], xa+1.5*d, apply(lam4, 2, qu)[2,], col=cl[4])
axis(2)
axis(1, at=1:6, labels=1988:1993)
legend(x=2.5, y=0.75, pch=rep(16,2), col=cl[1:2],
    legend=c(expression('P'[1]),expression('P'[2])), bty='n')
legend(x=3.5, y=0.75, pch=rep(16,2), col=cl[3:4], legend=c(expression('P'[3]), 'U'), bty='n')
par(op)


# Define an IPM that uses a uniform prior for the initial population size. In principle this is the same model as model2.txt, however, instead of using the categorical distribution in JAGS to create a discrete uniform prior, we use here the uniform distribution and round the simulated number to the nearest integer. This is computationally more efficient, and is expected to have no impact on the results.

# Write JAGS model file
cat(file = "model3.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.years-1)){
    f[t] ~ dunif(0, 1)
  }
  for (t in 1:n.years){
    s[t] ~ dunif(0, 1)
    h[t] ~ dunif(0, 1)
    r[t] ~ dunif(0, 1)
  }

  # Age-at-harvest data (state-space model)
  # Model for the initial population size: rounded uniform priors
  for (a in 1:n.age){
    n[a] ~ dunif(0.5, upper+0.4999)    # Ensures that 1 and upper are chosen with the same probability as values 2 to 999
    N[a,1] <- round(n[a])
  }

  # Process model over time: our model of population dynamics
  for (t in 1:(n.years-1)){
    N[1,t+1] ~ dpois((Ntot[t] - N[1,t]) * f[t])
    for (a in 2:n.age){
      N[a,t+1] ~ dbin((1-h[t]) * s[t], N[a-1,t])
    } #a
  } #t

  # Derived quantity: total year-specific population size
  for (t in 1:n.years){
    Ntot[t] <- sum(N[,t])
  }

  # Observation model
  for (t in 1:n.years){
    for (a in 1:n.age){
      C[a,t] ~ dbin(h[t] * r[t], N[a,t])
    } #a
  } #t

  # Hunter survey data (logistic regression model)
  for (t in 1:n.years){
    b[t] ~ dbin(r[t], a[t])
  }

  # Radio tracking data (multinomial model)
  for (t in 1:n.years){
    R[,t] ~ dmulti(prt[,t], total[t])
    prt[1,t] <- h[t]
    prt[2,t] <- (1-h[t]) * (1-s[t])
    prt[3,t] <- (1-h[t]) * s[t]
  }
}
")



# Fit models

# 1. Original data set; uniform prior U(1, 1000)
jags.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R, total=colSums(R),
    n.years=ncol(C), n.age=nrow(C), upper=1000))

# Initial values
# We need good initial values for the population size in the first year
s <- 0.8   # guess of annual survival
f <- 0.25  # guess of recruitment

# Create matrix population model
A <- matrix(0, ncol=17, nrow=17)
A[1,2:17] <- f
for (a in 2:17){
  A[a,a-1] <- s
}
# Compute stable age distribution (right eigenvector)
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])
# Population size in first age class of first year
r <- 0.5   # guess of reporting rate
h <- 0.1   # guess of hunting mortality
N1 <- elk$C[1,1] / (h * r)
# Compute age-specific population sizes in first year
n <- N1 * revec / revec[1]

inits <- function() {list(s=runif(jags.data$n.years, 0.8, 1), n=round(n),
    h=runif(jags.data$n.years, 0.05, 0.1))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
# ni <- 250000; nt <- 20; nb <- 50000; nc <- 3; na <- 5000
ni <- 25000; nt <- 2; nb <- 5000; nc <- 3; na <- 500  # ~~~ for testing

# Call JAGS from R (ART 4 min) and check convergence
out5 <- jags(jags.data, inits, parameters, "model3.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = na, parallel = TRUE)


# 2. Original data set; uniform prior U(1, 3000)
jags.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R, total=colSums(R),
    n.years=ncol(C), n.age=nrow(C), upper=3000))

# Initial values
inits <- function() {list(s=runif(jags.data$n.years, 0.8, 1), n=round(n),
    h=runif(jags.data$n.years, 0.05, 0.1))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
# ni <- 1600000; nt <- 150; nb <- 100000; nc <- 3; na <- 5000
ni <- 160000; nt <- 15; nb <- 10000; nc <- 3; na <- 500  # ~~~ for testing

# Call JAGS from R (ART 27 min) and check convergence
out6 <- jags(jags.data, inits, parameters, "model3.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = na, parallel = TRUE)


# 3. 20 times more radio tracking data; uniform prior U(1, 1000)
jags.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R*20, total=colSums(R*20),
    n.years=ncol(C), n.age=nrow(C), upper=1000))

# Initial values
inits <- function() {list(s=runif(jags.data$n.years, 0.8, 1), n=round(n),
    h=runif(jags.data$n.years, 0.05, 0.1))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
# ni <- 250000; nt <- 20; nb <- 50000; nc <- 3; na <- 5000
ni <- 25000; nt <- 2; nb <- 5000; nc <- 3; na <- 500  # ~~~ for testing

# Call JAGS from R (ART 4 min) and check convergence
out7 <- jags(jags.data, inits, parameters, "model3.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = na, parallel = TRUE)


# 4. 20 times more radio tracking data; uniform prior U(1, 3000)
jags.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R*20, total=colSums(R*20),
    n.years=ncol(C), n.age=nrow(C), upper=3000))

# Initial values
inits <- function() {list(s=runif(jags.data$n.years, 0.8, 1), n=round(n),
    h=runif(jags.data$n.years, 0.05, 0.1))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
# ni <- 250000; nt <- 20; nb <- 50000; nc <- 3; na <- 5000
ni <- 25000; nt <- 2; nb <- 5000; nc <- 3; na <- 500  # ~~~ for testing

# Call JAGS from R (ART 4 min) and check convergence
out8 <- jags(jags.data, inits, parameters, "model3.txt",
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = na, parallel = TRUE)


save(out1, out2, out3, out4, out5, out6, out7, out8, file="ElkResults.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~ Fig. 17.7 ~~~~
load('ElkResults.Rdata')
library(scales)
cl <- c(viridis_pal(option='E')(20)[c(18,8)], 'red')

op <- par("mfrow", "mar")
layout(matrix(1:16, 4, 4, byrow=TRUE), widths=c(1.1, 1, 1, 1), heights=c(1.1, 1, 1, 1), TRUE)

up <- 1000
upy <- 0.006

par(mar=c(3,3,2,1))
hist(out5$sims.list$N[,1,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
abline(h=1/up, col=cl[3])
text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out5$mean$N[1,1]))), pos=4)
text(x=0, y=0.8*upy, labels=paste('SD: ', round(out5$sd$N[1,1])), pos=4)
axis(1)
mtext('U(1,1000)', side=2, line=1)
mtext('1y old', side=3, line=0.5)

par(mar=c(3,2,2,1))
hist(out5$sims.list$N[,2,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
abline(h=1/up, col=cl[3])
text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out5$mean$N[2,1]))), pos=4)
text(x=0, y=0.8*upy, labels=paste('SD: ', round(out5$sd$N[2,1])), pos=4)
axis(1)
mtext('2y old', side=3, line=0.5)

par(mar=c(3,2,2,1))
hist(out5$sims.list$N[,3,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
abline(h=1/up, col=cl[3])
text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out5$mean$N[3,1]))), pos=4)
text(x=0, y=0.8*upy, labels=paste('SD: ', round(out5$sd$N[3,1])), pos=4)
axis(1)
mtext('3y old', side=3, line=0.5)

par(mar=c(3,2,2,1))
hist(out5$sims.list$N[,4,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
abline(h=1/up, col=cl[3])
text(x=600, y=0.9*upy, labels=bquote(bar(x) == .(round(out5$mean$N[4,1]))), pos=4)
text(x=600, y=0.8*upy, labels=paste('SD: ', round(out5$sd$N[4,1])), pos=4)
axis(1)
mtext('4y old', side=3, line=0.5)

par(mar=c(3,3,1,1))
up <- 3000
upy <- 0.0015
hist(out6$sims.list$N[,1,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
abline(h=1/up, col=cl[3])
text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out6$mean$N[1,1]))), pos=4)
text(x=0, y=0.8*upy, labels=paste('SD: ', round(out6$sd$N[1,1])), pos=4)
axis(1)
mtext('U(1,3000)', side=2, line=1)

par(mar=c(3,2,1,1))
hist(out6$sims.list$N[,2,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
abline(h=1/up, col=cl[3])
text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out6$mean$N[2,1]))), pos=4)
text(x=0, y=0.8*upy, labels=paste('SD: ', round(out6$sd$N[2,1])), pos=4)
axis(1)

par(mar=c(3,2,1,1))
hist(out6$sims.list$N[,3,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
abline(h=1/up, col=cl[3])
text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out6$mean$N[3,1]))), pos=4)
text(x=0, y=0.8*upy, labels=paste('SD: ', round(out6$sd$N[3,1])), pos=4)
axis(1)

par(mar=c(3,2,1,1))
hist(out6$sims.list$N[,4,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
abline(h=1/up, col=cl[3])
text(x=1700, y=0.9*upy, labels=bquote(bar(x) == .(round(out6$mean$N[4,1]))), pos=4)
text(x=1700, y=0.8*upy, labels=paste('SD: ', round(out6$sd$N[4,1])), pos=4)
axis(1)

par(mar=c(3,3,1,1))
up <- 1000
upy <- 0.008
hist(out7$sims.list$N[,1,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
abline(h=1/up, col=cl[3])
text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out7$mean$N[1,1]))), pos=4)
text(x=0, y=0.8*upy, labels=paste('SD: ', round(out7$sd$N[1,1])), pos=4)
axis(1)
mtext('U(1,1000)', side=2, line=1)

par(mar=c(3,2,1,1))
hist(out7$sims.list$N[,2,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
abline(h=1/up, col=cl[3])
text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out7$mean$N[2,1]))), pos=4)
text(x=0, y=0.8*upy, labels=paste('SD: ', round(out7$sd$N[2,1])), pos=4)
axis(1)

par(mar=c(3,2,1,1))
hist(out7$sims.list$N[,3,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
abline(h=1/up, col=cl[3])
text(x=500, y=0.9*upy, labels=bquote(bar(x) == .(round(out7$mean$N[3,1]))), pos=4)
text(x=500, y=0.8*upy, labels=paste('SD: ', round(out7$sd$N[3,1])), pos=4)
axis(1)

par(mar=c(3,2,1,1))
hist(out7$sims.list$N[,4,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
abline(h=1/up, col=cl[3])
text(x=500, y=0.9*upy, labels=bquote(bar(x) == .(round(out7$mean$N[4,1]))), pos=4)
text(x=500, y=0.8*upy, labels=paste('SD: ', round(out7$sd$N[4,1])), pos=4)
axis(1)

par(mar=c(3,3,1,1))
up <- 3000
upy <- 0.006
hist(out8$sims.list$N[,1,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
abline(h=1/up, col=cl[3])
text(x=1000, y=0.9*upy, labels=bquote(bar(x) == .(round(out8$mean$N[1,1]))), pos=4)
text(x=1000, y=0.8*upy, labels=paste('SD: ', round(out8$sd$N[1,1])), pos=4)
axis(1)
mtext('U(1,3000)', side=2, line=1)

par(mar=c(3,2,1,1))
hist(out8$sims.list$N[,2,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
abline(h=1/up, col=cl[3])
text(x=1000, y=0.9*upy, labels=bquote(bar(x) == .(round(out8$mean$N[2,1]))), pos=4)
text(x=1000, y=0.8*upy, labels=paste('SD: ', round(out8$sd$N[2,1])), pos=4)
axis(1)

par(mar=c(3,2,1,1))
hist(out8$sims.list$N[,3,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
abline(h=1/up, col=cl[3])
text(x=1500, y=0.9*upy, labels=bquote(bar(x) == .(round(out8$mean$N[3,1]))), pos=4)
text(x=1500, y=0.8*upy, labels=paste('SD: ', round(out8$sd$N[3,1])), pos=4)
axis(1)

par(mar=c(3,2,1,1))
hist(out8$sims.list$N[,4,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
    ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
abline(h=1/up, col=cl[3])
text(x=1500, y=0.9*upy, labels=bquote(bar(x) == .(round(out8$mean$N[4,1]))), pos=4)
text(x=1500, y=0.8*upy, labels=paste('SD: ', round(out8$sd$N[4,1])), pos=4)
axis(1)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 17.8 Specification of the survival process with hazard rates
# ============================================================

# ~~~~ Recreate the priors used for the first IPM (lines 36-56 above ~~~~
s <- 0.8                    # Guess of annual survival
f <- 0.25                   # Guess of recruitment

# Create matrix population model
A <- matrix(0, ncol=17, nrow=17)
A[1,2:17] <- f
for (a in 2:17){
  A[a,a-1] <- s
}

# Compute right eigenvector
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])

# Population size in first age class of first year
r <- 0.5                    # Guess of reporting rate
h <- 0.1                    # Guess of hunting mortality
N1 <- elk$C[1,1] / (h * r)

# Compute age-dependent population size in first year
n <- N1 * revec / revec[1]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Bundle data and produce data overview
jags.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R, total=colSums(R), n.years=ncol(C),
    n.age=nrow(C), n=n))

# Write JAGS model file
cat(file = "model3.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.years-1)){
    f[t] ~ dunif(0, 1)
  }
  for (t in 1:n.years){
    # Priors for hazard rates
    mh[t] ~ dgamma(0.1, 0.1)                              # Hunting mortality hazard rate
    mo[t] ~ dgamma(0.1, 0.1)                              # Background mortality hazard rate
    # Calculate probabilities from hazard rates
    s[t] <- exp(-(mh[t] + mo[t]))                         # Overall survival
    h[t] <- (1 - s[t]) * (mh[t] / (mh[t] + mo[t]))        # Hunting mortality
    o[t] <- (1 - s[t]) * (mo[t] / (mh[t] + mo[t]))        # Background mortality
    # Prior for reporting rate
    r[t] ~ dunif(0, 1)
  }

  # Age-at-harvest data (state-space model)
  # Model for the initial population size: Poisson priors
  for (a in 1:n.age){
    N[a,1] ~ dpois(n[a])
  }

  # Process model over time: our model of population dynamics
  for (t in 1:(n.years-1)){
    N[1,t+1] ~ dpois((Ntot[t] - N[1,t]) * f[t])
    for (a in 2:n.age){
      N[a,t+1] ~ dbin(s[t], N[a-1,t])
    } #a
  } #t

  # Derived quantity: total year-specific population size
  for (t in 1:n.years){
    Ntot[t] <- sum(N[,t])
  }

  # Observation model
  for (t in 1:n.years){
    for (a in 1:n.age){
      C[a,t] ~ dbin(h[t] * r[t], N[a,t])
    } #a
  } #t

  # Hunter survey data (logistic regression model)
  for (t in 1:n.years){
    b[t] ~ dbin(r[t], a[t])
  }

  # Radio tracking data (multinomial model)
  for (t in 1:n.years){
    R[,t] ~ dmulti(prt[,t], total[t])
    prt[1,t] <- h[t]
    prt[2,t] <- o[t]
    prt[3,t] <- s[t]
  }
}
")

# Initial values
inits <- function() {list(mh=runif(jags.data$n.years, 0.01, 0.1), mo=runif(jags.data$n.years,
    0.01, 0.1))}

# Parameters monitored
parameters <- c("mh", "mo", "h", "s", "o", "f", "r", "Ntot", "N")

# MCMC settings
ni <- 450000; nb <- 50000; nc <- 3; nt <- 40; na <- 5000

# Call JAGS from R (ART 9 min) and check convergence
out9 <- jags(jags.data, inits, parameters, "model3.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out9)
print(out9, 3)

#              mean      sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# mh[1]       0.129   0.015    0.102    0.128    0.161    FALSE 1 1.001  8101
# mh[2]       0.152   0.020    0.118    0.150    0.197    FALSE 1 1.003  1669
# mh[3]       0.083   0.011    0.065    0.082    0.110    FALSE 1 1.004  1791
# mh[4]       0.078   0.012    0.058    0.076    0.106    FALSE 1 1.004  1279
# mh[5]       0.075   0.012    0.056    0.073    0.104    FALSE 1 1.005  1123
# mh[6]       0.174   0.038    0.118    0.169    0.265    FALSE 1 1.004  1287
# mo[1]       0.112   0.071    0.016    0.098    0.283    FALSE 1 1.004  1132
# mo[2]       0.039   0.036    0.001    0.029    0.133    FALSE 1 1.000 30000
# mo[3]       0.038   0.036    0.001    0.028    0.134    FALSE 1 1.001 30000
# mo[4]       0.004   0.010    0.000    0.000    0.036    FALSE 1 1.192    64
# mo[5]       0.173   0.083    0.049    0.160    0.371    FALSE 1 1.001  6370
# mo[6]       0.189   0.094    0.052    0.175    0.412    FALSE 1 1.000 15474
# h[1]        0.114   0.012    0.092    0.114    0.139    FALSE 1 1.000 30000
# h[2]        0.138   0.017    0.109    0.137    0.176    FALSE 1 1.002  1774
# h[3]        0.079   0.010    0.062    0.077    0.102    FALSE 1 1.004  1949
# h[4]        0.075   0.011    0.056    0.073    0.101    FALSE 1 1.004  1174
# h[5]        0.066   0.010    0.050    0.065    0.090    FALSE 1 1.005  1314
# h[6]        0.146   0.028    0.103    0.142    0.212    FALSE 1 1.004  1159
# s[1]        0.788   0.058    0.655    0.797    0.874    FALSE 1 1.003  1236
# s[2]        0.827   0.035    0.742    0.834    0.877    FALSE 1 1.001  5186
# s[3]        0.886   0.034    0.801    0.894    0.928    FALSE 1 1.001  7629
# s[4]        0.922   0.015    0.885    0.924    0.943    FALSE 1 1.041   166
# s[5]        0.783   0.066    0.637    0.790    0.890    FALSE 1 1.001  4348
# s[6]        0.699   0.072    0.540    0.705    0.819    FALSE 1 1.000 27082
# o[1]        0.097   0.057    0.015    0.088    0.231    FALSE 1 1.003  1230
# o[2]        0.035   0.031    0.001    0.026    0.116    FALSE 1 1.000 30000
# o[3]        0.035   0.032    0.001    0.026    0.120    FALSE 1 1.001 30000
# o[4]        0.003   0.010    0.000    0.000    0.034    FALSE 1 1.188    64
# o[5]        0.150   0.066    0.046    0.143    0.299    FALSE 1 1.000  7247
# o[6]        0.155   0.068    0.047    0.147    0.309    FALSE 1 1.000 11278
# f[1]        0.271   0.036    0.203    0.270    0.343    FALSE 1 1.002  1744
# f[2]        0.312   0.038    0.242    0.311    0.391    FALSE 1 1.001  2894
# f[3]        0.327   0.043    0.250    0.325    0.418    FALSE 1 1.001  1956
# f[4]        0.219   0.034    0.157    0.217    0.292    FALSE 1 1.000  5287
# f[5]        0.146   0.035    0.087    0.143    0.223    FALSE 1 1.000 29143
# r[1]        0.522   0.030    0.464    0.522    0.581    FALSE 1 1.000  6106
# r[2]        0.531   0.029    0.473    0.530    0.587    FALSE 1 1.000 30000
# r[3]        0.995   0.005    0.983    0.997    1.000    FALSE 1 1.000 30000
# r[4]        0.755   0.023    0.710    0.756    0.798    FALSE 1 1.000 30000
# r[5]        0.995   0.005    0.981    0.997    1.000    FALSE 1 1.000 22544
# r[6]        0.481   0.027    0.428    0.481    0.535    FALSE 1 1.000 30000
# Ntot[1]  2097.226  44.984 2010.000 2097.000 2186.000    FALSE 1 1.000 30000
# Ntot[2]  2093.870 165.500 1727.000 2111.000 2371.000    FALSE 1 1.003  1318
# Ntot[3]  2236.494 225.383 1756.000 2252.000 2637.025    FALSE 1 1.004  1308
# Ntot[4]  2537.690 295.944 1917.000 2553.000 3080.000    FALSE 1 1.004  1048
# Ntot[5]  2761.377 352.628 2031.975 2776.000 3406.025    FALSE 1 1.003  1260
# Ntot[6]  2491.857 398.151 1714.000 2491.000 3268.000    FALSE 1 1.002  1298
# N[1,1]    428.932  18.692  392.000  429.000  466.000    FALSE 1 1.000 30000
# N[2,1]    355.061  16.618  323.000  355.000  388.000    FALSE 1 1.000 30000
# [...output truncated...]
# N[16,6]    19.665   4.792   11.000   19.000   30.000    FALSE 1 1.000  5106
# N[17,6]    16.137   4.328    8.000   16.000   25.000    FALSE 1 1.000  5597
