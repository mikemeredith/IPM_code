# Schaub & Kéry (2022) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------

# Run time approx. 8 mins

library(IPMbook) ; library(jagsUI)

# 4.6 Brief introduction to spatial capture-recapture (SCR) modeling
# ==================================================================

# ~~~ This does not work with current versions of R ~~~
# library(scrbook)                      # From sites.google.com/site/spatialcapturerecapture/
# ?simSCR0                              # Open help page

# str(dat <- simSCR0(N=100, K=20, discard0=TRUE, array3d=FALSE, rnd=2013))

N <- 100                                # Population size: number of individuals living in state space
K <- 20                                 # Number of trap nights
rnd <- 2013                             # Random number seed

traplocs <- cbind(sort(rep(1:5, 5)), rep(1:5, 5)) # Set up trapping grid
J <- nrow(traplocs)                               # Number of traps (ntraps)

# Set up state space
buffer <- 2                             # Buffer size = 2 units
Xl <- min(traplocs[, 1] - buffer)       # X lower
Xu <- max(traplocs[, 1] + buffer)       # X upper
Yl <- min(traplocs[, 2] - buffer)       # Y lower
Yu <- max(traplocs[, 2] + buffer)       # Y upper

set.seed(rnd)                           # Set seed
sx <- runif(N, Xl, Xu)
sy <- runif(N, Yl, Yu)
S <- cbind(sx, sy)
head(S)                                 # Look at coordinates of activity centers (ACs) - not shown

library(AHMbook)                        # In case scrbook does not work
D <- AHMbook::e2dist(S, traplocs)       # Distances between AC and traps
str(D)                                  # Shows a 100 x 25 matrix

# Plot realization of spatial point process: ACs and trapping grid
# and variability in the distance between AC and traps (Fig. 4.15)
plot(traplocs, xlim=c(Xl, Xu), ylim=c(Yl, Yu), pch=16, axes=FALSE, xlab='x-coordinate',
    ylab='y-coordinate', cex=1)
points(S, col="red", pch=1, cex=0.8)
axis(2); axis(1)

# ~~~~ histogram of distances ~~~~
hist(D, breaks=50, col="grey", main="", xlab='Distance between trap and AC',
    axes=FALSE, border=FALSE)
axis(2); axis(1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 4.16 left: relationship between p and scale parameter sigma ~~~~
alpha0 <- -2.5
sigma <- 0.5
op <- par(mfrow=c(1, 2), mar=c(5, 5, 4, 3), cex.lab=1.5, cex.axis=1.5)
SIGMA <- round(seq(0.1, 10, , 8), 2)    # Choose 8 values of sigma between 0.1 and 10
xd <- seq(0, 8, , 1000)
plot(xd, plogis(alpha0) * exp(-(1/(2 * SIGMA[1]^2)) * xd^2),
    xlab="Distance from trap to activity center (AC)",
    ylab="Detection probability (p)", type="l", col="grey",
    frame.plot=FALSE, lwd=3)
for(k in 2:8){
  lines(xd, plogis(alpha0) * exp(-(1/(2 * SIGMA[k]^2)) * xd^2), type="l", col="grey", lwd=3)
}
lines(xd, plogis(alpha0) * exp(-(1/(2 * sigma^2)) * xd^2), type="l", col="black", lwd=3)
text(1.2, 0.021, labels=sigma, cex=1)
text(0.45, 0.02, labels=round(SIGMA[1],1), cex=1, col="grey")
text(2.6, 0.025, labels=round(SIGMA[2],1), cex=1, col="grey")
text(4.3, 0.032, labels=round(SIGMA[3],1), cex=1, col="grey")
text(5.6, 0.038, labels=round(SIGMA[4],1), cex=1, col="grey")
text(6.5, 0.044, labels=round(SIGMA[5],1), cex=1, col="grey")
text(7.2, 0.05, labels=round(SIGMA[6],1), cex=1, col="grey")
text(7.6, 0.054, labels=round(SIGMA[7],1), cex=1, col="grey")
text(7.8, 0.059, labels=round(SIGMA[8],1), cex=1, col="grey")
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define parameters of detection probability
alpha0 <- -2.5                          # Logit(baseline detection probability right on a trap)
plogis(alpha0)                          # p0: detection probability for AC right ON a trap: 0.0758
sigma <- 0.5                            # Scale parameter of half-normal kernel
alpha1 <- 1 / (2 * sigma^2)             # Coefficient on distance

# Compute probability of detection for each combination of animal and trap
probcap <- plogis(alpha0) * exp(-alpha1 * D^2)
str(probcap)                            # This is a 100 x 25 matrix

# ~~~~ Fig. 4.16 right: Frequency distribution of capture probs ~~~~
hist(probcap, col="dodgerblue", main="", xlim=c(0.001,0.09), freq=FALSE,
    breaks=40, xlab='Detection probability')
abline(v=plogis(alpha0), col="black", lwd=2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ What is p at distance sigma, for sigma = 0.5 ? ~~~~
plogis(alpha0) * exp(-alpha1 * sigma^2)
# [1] 0.04601031
exp(-alpha1 * sigma^2)                  # Drop to about 60% of max(p), which is p0
# [1] 0.6065307
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

y3d <- array(NA, dim=c(N, J, K))                  # Dimension is individual x trap x night
for (i in 1:nrow(y3d)){                           # Loop over individuals
  for (j in 1:J){                                 # Loop over traps
    y3d[i,j,] <- rbinom(K, 1, probcap[i,j])       # Bernoulli trials
  } #j
} #i

# ~~~~ Using scrbook::simSCR0 ~~~~
# As said before, you can greatly abbreviate the generation of such data sets by using function simSCR0 (though the choice of seed will be screwed):

if(require(scrbook, quietly = TRUE)) {
  str(dat <- simSCR0(N=100, K=20, discard0=TRUE, array3d=FALSE, rnd=2013))
}
# List of 9
# $ Y       : int [1:42, 1:25] 0 0 0 0 0 0 0 0 0 0 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:42] "1" "2" "3" "4" ...
# .. ..$ : chr [1:25] "trap1" "trap2" "trap3" "trap4" ...
# $ traplocs: int [1:25, 1:2] 1 1 1 1 1 2 2 2 2 2 ...
# $ xlim    : num [1:2] -1 7
# $ ylim    : num [1:2] -1 7
# $ N       : num 100
# $ alpha0  : num -2.5
# $ alpha1  : num 2
# $ sigma   : num 0.5
# $ K       : num 20

# Setting discard0 = FALSE would return the complete detection matrix including for the animals never seen, while setting array3d = TRUE will produce data in the Bernoulli format, rather than as a binomial count aggregating over trap nights. The detection function parameters are currently hard-wired, but it would be easy to adapt the function such that they could be varied and included as function arguments. (This would actually be a great exercise for you. Writing functions in R is essential for you as a data analyst and extending an existing function to include additional functionality is what you will fairly often do.)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

y2d <- apply(y3d, c(1, 2), sum)         # Binomial detection frequencies
str(y2d)

table(y2d)
# y2d
#    0    1    2    3    4    5
# 2441   45   10    2    1    1

table(apply(y2d, 1, sum))
#  0  1  2  3  4  5  6
# 65 10 15  4  3  2  1

yy <- y2d
yy[yy > 1] <- 1
table(apply(yy, 1, sum))
#  0  1  2  3  4
# 65 17 13  4  1

totalcaps <- apply(y2d, 1, sum)
nobs <- sum(totalcaps > 0)
y <- y2d[totalcaps > 0, ]
dimnames(y) <- list(1:nrow(y), paste("trap", 1:ncol(y), sep=""))
str(y)
# int [1:35, 1:25] 0 0 0 0 1 0 0 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:35] "1" "2" "3" "4" ...
# ..$ : chr [1:25] "trap1" "trap2" "trap3" "trap4" ...

Sobs <- S[totalcaps > 0,]               # Activity center coordinates of detected animals
# op <- par(ask=TRUE)
op <- par(ask=dev.interactive(orNone=TRUE))     # ~~~ for testing
for (i in 1:nrow(y)){                           # Loop over individuals
  plot(traplocs, xlim=c(Xl, Xu), ylim=c(Yl, Yu), pch=16, xlab='x-coordinate', ylab='y-coordinate',
      frame=FALSE)
  points(S[totalcaps == 0,], col='red', pch=1)  # Never captured
  points(S[totalcaps > 0,], col='red', pch=16)  # Captured >0 times
  for (k in 1:i){
    tmp <- which(y[k,] > 0)                     # Detected by which trap
    for (j in 1:length(tmp)){
      segments(Sobs[k,1], Sobs[k,2], traplocs[tmp,1], traplocs[tmp,2], col='gray', lty=1)
    } #j
  } #k
  for (j in 1:length(tmp)){
    segments(Sobs[i,1], Sobs[i,2], traplocs[tmp,1], traplocs[tmp,2], col='red', lwd=2)
  } #j
  title(main=paste('Lynx', i, 'caught in: ', paste(names(tmp), collapse=", ")))
} #i
par(op)

nind <- nrow(y)                                         # Observed number of individuals
J <- nrow(traplocs)                                     # Number of traps
xlim <- c(Xl, Xu)                                       # x limits of state-space
ylim <- c(Yl, Yu)                                       # y limits of state-space
A <- (max(xlim) - min(xlim)) * (max(ylim) - min(ylim))  # State-space area

# Data augmentation
M <- 200                                                # Allow for up to M-35=165 undetected lynx
y <- rbind(y, matrix(0, nrow=M - nind, ncol=ncol(y)))

# Bundle data
jags.data <- list(y=y, traplocs=traplocs, K=K, M=M, J=J, xlim=xlim, ylim=ylim, A=A)
str(jags.data)
# List of 8
# $ y       : num [1:200, 1:25] 0 0 0 0 1 0 0 0 0 0 ...
# $ traplocs: int [1:25, 1:2] 1 1 1 1 1 2 2 2 2 2 ...
# $ K       : num 20
# $ M       : num 200
# $ J       : int 25
# $ xlim    : num [1:2] -1 7
# $ ylim    : num [1:2] -1 7
# $ A       : num 64

# Write JAGS model file
cat(file="model23.txt", "
model {
  # Priors and linear models
  p0 ~ dunif(0, 1)                      # Baseline detection probability
  alpha ~ dnorm(0, 0.1)                 # 1 / (2*sigma^2) (coef of distance)
  sigma <- sqrt(1 / (2 * alpha))        # Half-normal scale
  psi ~ dunif(0, 1)                     # Data augmentation parameter

  # Likelihood
  for (i in 1:M){                       # Loop over all M=200 individuals
    z[i] ~ dbern(psi)                   # Data augmentation variable
    # State model: point process model
    s[i,1] ~ dunif(xlim[1], xlim[2])    # x-coord of activity centers
    s[i,2] ~ dunif(ylim[1], ylim[2])    # y coord of activity centers
    # Observation model: p ~ distance between trap and estimated AC
    for (j in 1:J){                     # Loop over all traps
      d2[i,j] <- pow(s[i,1] - traplocs[j,1], 2) + pow(s[i,2] - traplocs[j,2], 2) # Distance squared
      p[i,j] <- z[i] * p0 * exp(-alpha * d2[i,j]) # Detection prob.
      y[i,j] ~ dbin(p[i,j],K)           # Observed detection frequencies
    } #j
  } #i

  # Derived quantities
  N <- sum(z[])                         # Population size in state-space
  D <- N / A                            # Density over state-space
}
")

# Function to generate starting values
sst <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
for (i in 1:nind){
  if (sum(y[i, ]) == 0)
    next
  sst[i, 1] <- mean(traplocs[y[i, ] > 0, 1])
  sst[i, 2] <- mean(traplocs[y[i, ] > 0, 2])
}
zst <- c(rep(1, nind), rep(0, M - nind))
inits <- function() {list(p0=runif(1), alpha=runif(1, 1, 2), s=sst, z=zst)}

# Parameters to save
parameters <- c("alpha", "sigma", "N", "D", "s", "z")

# MCMC settings
ni <- 15000; nb <- 5000; nc <- 3; nt <- 5; na <- 1000

# Call JAGS (ART 7 min), check convergence and summarize posteriors
out26 <- jags(jags.data, inits, parameters, "model23.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out26) # Not shown
print(out26, 3)
#          mean     sd   2.5%    50%   97.5% overlap0     f  Rhat n.eff
# alpha   1.878  0.308  1.323  1.860   2.533    FALSE 1.000 1.002  1022
# sigma   0.521  0.043  0.444  0.519   0.615    FALSE 1.000 1.002  1102
# N      78.992 11.405 59.000 78.000 104.000    FALSE 1.000 1.001  1236
# D       1.234  0.178  0.922  1.219   1.625    FALSE 1.000 1.001  1236
# s[1,1]  3.006  0.377  2.258  3.007   3.740    FALSE 1.000 1.000  6000
# s[2,1]  4.720  0.331  4.091  4.706   5.406    FALSE 1.000 1.000  4089
# s[3,1]  5.226  0.408  4.416  5.242   5.972    FALSE 1.000 1.001  2284
# [... output truncated ...]


# ~~~~ Plot posteriors of population size and density parameters, figure 4.18 ~~~~
op <- par(mfrow=c(1, 2), mar=c(5,5,4,3), cex.axis=1.25, cex.lab=1.25, las=1)
hist(out26$sims.list$N, breaks=40, col="grey",
    xlab=expression(paste("Population size (", italic(N), ")")),
    xlim=range(c(20, nobs, max(out26$sims.list$N))),
    freq=FALSE, main="", border=NA, lwd=2)
abline(v=nobs, col="black", lwd=3)
abline(v=out26$mean$N, col="blue", lwd=3)
abline(v=N, col="red", lwd=3)

hist(out26$sims.list$D, breaks=40, col="grey",
    xlab=expression(paste("Population density (", italic(D), ")")),
    xlim=c(0.5,2.5), freq=FALSE, main="", border=NA, ylab=NA, lwd=2)
abline(v=nobs/16, col="black", lwd=3)
abline(v=out26$mean$D, col="blue", lwd=3)
abline(v=100/A, col="red", lwd=3)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for figure 4.19 ~~~~
op <- par(mfrow=c(1, 2), mar= c(5, 5, 4, 3), cex.lab=1.5, cex.axis=1.5)
# Fig. 4.19 (left): ACs of the 35 animals detected at least once
plot(traplocs, xlim=c(Xl, Xu), ylim=c(Yl, Yu), pch=16, cex=1,
    frame=FALSE, xlab='x-coordinate', ylab='y-coordinate', asp=1)
points(S[totalcaps == 0,], col="red", pch=1, cex=1) # Those never captured
points(S[totalcaps > 0,], col="red", pch=16, cex=1) # Those captured at least once
for (i in 1:nobs){
  points(out26$sims.list$s[,i,], col=i, pch='.')    # Estimated activity centers
}
points(traplocs, pch=16, cex.main=1)

# Fig. 4.19 (right): ACs of the 65 animals never detected
plot(traplocs, xlim=c(Xl, Xu), ylim=c(Yl, Yu), pch=16, cex=1,
    frame=FALSE, xlab='x-coordinate', ylab='y-coordinate', asp=1)
points(S[totalcaps == 0,], col="red", pch=1, cex=1) # Those never captured
points(S[totalcaps > 0,], col="red", pch=16, cex=1) # Those captured at least once
for (i in (nobs+1):M){
  real.ind <- out26$sims.list$z[,i] == 1
  points(out26$sims.list$s[real.ind,i,], col=i, pch='.')
}
points(traplocs, pch=16, cex.main=1)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Produce a density map based on the estimated ACs (Fig. 4.20 left)
op <- par(mfrow=c(1, 2))
cl <- terrain.colors(20, rev=TRUE)
# Extract posterior samples for s and z
z <- out26$sims.list$z
s <- out26$sims.list$s
s[z == 0] <- NA                         # 'Remove' ghost individuals
niter <- nrow(z)
# Choose a grid
xg <- seq(xlim[1], xlim[2], , 50)
yg <- seq(ylim[1], ylim[2], , 50)
# Associate each draw with its pixel
Lx1 <- cut(s[, , 1], breaks=xg, include.lowest=TRUE)
Ly1 <- cut(s[, , 2], breaks=yg, include.lowest=TRUE)
# Tally up the number of samples per pixel: a density table
D.hat1 <- table(Lx1, Ly1) / niter
sum(D.hat1)                             # Check, should be the same as the posterior mean of N

# Make a plot of that density table
image(xg, yg, D.hat1, col=cl, xlab='x-coordinate', ylab='y-coordinate', asp=1,
    main='Density based on home-range centers')
points(Sobs, col='red', pch=16)
points(traplocs, pch=16)
points(S, col='red')

# Produce a density map based on the estimated locations (Fig. 4.20 right)
# Collect up the stuff we need
AC <- out26$sims.list$s
sigma <- out26$sims.list$sigma
nind <- ncol(z)
# Do the posterior predictions of location
u <- array(NA, dim=dim(AC))
set.seed(1)
for (i in 1:niter){
  # Generate coordinates from normal distribution
  x <- rnorm(nind, AC[i,,1], sigma[i])
  y <- rnorm(nind, AC[i,,2], sigma[i])
  # Animals outside state space get back into the fold
  x <- ifelse(x > Xu, 2*Xu - x, x)
  x <- ifelse(x < Xl, 2*Xl - x, x)
  y <- ifelse(y > Xu, 2*Yu - y, y)
  y <- ifelse(y < Xl, 2*Yl - y, y)
  u[i,,] <- cbind(x, y)
}
u[z == 0] <- NA                         # 'Remove' ghost individuals
# Associate each draw with its pixel
Lx2 <- cut(u[,,1], breaks=xg, include.lowest=TRUE)
Ly2 <- cut(u[,,2], breaks=yg, include.lowest=TRUE)
# Tally up the number of samples per pixel: a density table
D.hat2 <- table(Lx2, Ly2) / niter
sum(D.hat2)                             # Check, should be the same as the posterior mean of N
# Make a plot of that density table
image(xg, yg, D.hat2, col=cl, xlab='x-coordinate', ylab='y-coordinate', asp=1,
    main='Density based on animal locations')
points(Sobs, col='red', pch=16)
points(traplocs, pch=16)
points(S, col='red')
par(op)
