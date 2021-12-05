# Schaub & Kéry (2022) Integrated Population Models
# Chapter 21 : Black bear
# -----------------------

# Run time for this test script approx 1.6 hrs, full run 15 hrs

# 21.4 Component data likelihoods
# ===============================

library(IPMbook); library(jagsUI)
data(bear)
str(bear)
# List of 3
# $ scr       : int [1:78, 1:128, 1:8, 1:3] 0 0 0 0 0 0 0 0 0 0 ...
# $ occ       : num [1:128, 1:8, 1:3] 0 0 0 0 0 NA 1 1 0 0 ...
# $ trap.coord: num [1:128, 1:2] 618 619 618 620 621 ...


# 21.4.1 Spatial capture-recapture data
# -------------------------------------

library(scales)
plot(bear$trap.coord, xlab='x coordinate', ylab='y coordinate', frame=FALSE, las=1, pch=16,
    col='#002A64FF', asp=1)
# Define limits of the state-space and add it to the plot
xlims <- c(615, 645)
ylims <- c(3380, 3420)
rect(xlims[1], ylims[1], xlims[2], ylims[2], col=alpha('grey', 0.3), border=NA)

table(bear$scr)                                   # Total number of detections
#      0   1
# 174640 626

# Replace the missing values by zero
bear.sum <- bear$scr
bear.sum[is.na(bear.sum)] <- 0

# Produce table
bearXtrap <- apply(bear.sum, c(1, 2), sum)

table(bearXtrap)
#    0    1    2    3    4    5    6    7    8    9   10
# 9638  209   82   23   10    6    8    3    2    1    2

table(apply(bearXtrap, 1, sum))
#  1  2  3  4  5  6  7  8  9 11 13 14 18 20 21 23 25 27 38 91
# 10  7 11  7  9  5  7  3  6  1  1  1  2  2  1  1  1  1  1  1


# Loop over 78 individuals and map their locations (individually)
# combined over all 3 years (which represent a 5-year duration)
# op <- par(ask=TRUE)                               # Press Esc to interrupt
op <- par(ask=dev.interactive(orNone=TRUE))         # Only ask if plotting on-screen
for (i in 1:nrow(bearXtrap)){
  cat(paste("\n\n*** Plot bear number", i, "***\n\n"))
  lucky.traps <- which(bearXtrap[i,] > 0)
  plot(bear$trap.coord, xlab='x coordinate', ylab='y coordinate', frame=FALSE, pch=1,
      col=rgb(0,0,0,1), cex=1.5, asp=1, main=paste('Bear number', i))
  points((bear$trap.coord)[lucky.traps,1], (bear$trap.coord)[lucky.traps,2], pch=16, col='red',
      cex=1.5)
}
par(op)


# 21.4.2 Occupancy data (no code)

# 21.5 The integrated population model
# ====================================

# Do data augmentation of the SCR data
M <- 200
y <- array(0, dim=c(M, dim(bear$scr)[2:4]))
y[1:dim(bear$scr)[1],,,] <- bear$scr

# Bundle data and produce data overview
jags.data <- list(y=y, o=bear$occ, T=6, K=dim(y)[3], J=dim(y)[2], M=M, x=bear$trap.coord, xlims=xlims,
    ylims=ylims)
str(jags.data)
# List of 9
# $ y    : num [1:200, 1:128, 1:8, 1:3] 0 0 0 0 0 0 0 0 0 0 ...
# $ o    : num [1:128, 1:8, 1:3] 0 0 0 0 0 NA 1 1 0 0 ...
# $ T    : num 6
# $ K    : int 8
# $ J    : int 128
# $ M    : num 200
# $ x    : num [1:128, 1:2] 618 619 618 620 621 ...
# $ xlims: num [1:2] 615 645
# $ ylims: num [1:2] 3380 3420


# Write JAGS model file
cat(file="model1.txt", "
model {
  # Priors
  psi ~ dunif(0, 1)                               # M*psi = E[N(1)], for data augmentation
  phi ~ dunif(0, 1)                               # Survival
  gamma ~ dunif(0, 3)                             # Per-capita recruitment
  p0 ~ dunif(0, 1)                                # Baseline capture probability
  sigma ~ dunif(0, 50)                            # Scale parameter of detection function

  # Derived parameters
  for (t in 1:T){
    N[t] <- sum(z[,t])                            # Population size at time t
    EB[t] <- N[t] * gamma                         # Expected number of recruits
    A[t] <- max(M - sum(a[,t]), 0.001)            # Bears available to be recruited
    b[t] <- min(EB[t] / A[t], 0.999)              # Probability of being recruited
  }

  # Population model (individual-based model)
  # Initial state
  for (i in 1:M){
    z[i,1] ~ dbern(psi)                           # Is a member of M alive?
    a[i,1] <- z[i,1]                              # Recruited yet?
    s[i,1] ~ dunif(xlims[1], xlims[2])            # Activity center, homogeneous
    s[i,2] ~ dunif(ylims[1], ylims[2])

    # Dynamics over time
    for (t in 2:T){
      z[i,t] ~ dbern(z[i,t-1] * phi + (1 - a[i,t-1]) * b[t-1])
      a[i,t] <- max(z[i,1:t])
    } #t
  } #i

  # Observation models
  # Spatial capture-recapture data
  for (i in 1:M){
    for (j in 1:J){
      d2[i,j] <- (s[i,1] - x[j,1])^2 + (s[i,2] - x[j,2])^2  # Distance squared
      p[i,j] <- p0 * exp(-d2[i,j] / (2 * sigma^2))
    } #j
    for (j in 1:J){                               # Trap
      for (k in 1:K){                             # Secondary occasions
        # The years with SCR data
        y[i,j,k,1] ~ dbern(p[i,j] * z[i,1])
        y[i,j,k,2] ~ dbern(p[i,j] * z[i,3])
        y[i,j,k,3] ~ dbern(p[i,j] * z[i,5])
      } #k
    } #j
    zi[i] <- (sum(z[i,]) > 0)                     # Was this bear ever alive?
  } #i

  # Occupancy data
  for (j in 1:J){
    for (k in 1:K){
      o[j,k,1] ~ dbern(1-prod(1-p[,j] * z[,2]))
      o[j,k,2] ~ dbern(1-prod(1-p[,j] * z[,4]))
      o[j,k,3] ~ dbern(1-prod(1-p[,j] * z[,6]))
    } #k
  } #j
  N.ever <- sum(zi[])                             # Bears ever alive – superpopulation size
}
")

# Initial values
zin <- array(0, dim=c(dim(y)[1],6))
zin[1:150,] <- 1
inits <- function() {list(z=zin, phi=runif(1, 0.7, 0.9), p0=runif(1, 0.01, 0.1),
    psi=runif(1, 0.2, 0.4), sigma=runif(1, 3, 4), gamma=runif(1, 0.2, 0.3))}

# Parameters monitored
parameters <- c("psi", "phi", "gamma", "p0", "sigma", "N", "EB", "b", "A", "s", "zi", "z", "N.ever")

# MCMC settings
# ni <- 5000; nb <- 2000; nc <- 3; nt <- 3; na <- 2000  # 14 hrs
ni <- 500; nb <- 200; nc <- 3; nt <- 1; na <- 200  # ~~~ for testing, 90 mins

# Call JAGS from R (ART 1400 min) and check convergence
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out1)

# 21.6 Results
# ============

print(out1, 3)

#              mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# psi         0.243  0.037    0.175    0.241    0.320    FALSE 1 1.001  3000
# phi         0.807  0.035    0.736    0.808    0.871    FALSE 1 1.001  1485
# gamma       0.235  0.037    0.166    0.233    0.312    FALSE 1 1.001  3000
# p0          0.115  0.007    0.101    0.115    0.129    FALSE 1 1.004   491
# sigma       3.205  0.093    3.042    3.198    3.404    FALSE 1 1.018   129
# N[1]       48.017  4.074   41.000   48.000   57.000    FALSE 1 1.003   665
# N[2]       59.722  4.278   52.000   60.000   69.000    FALSE 1 1.005   357
# N[3]       71.046  3.386   65.000   71.000   78.000    FALSE 1 1.008   248
# N[4]       64.342  4.704   56.000   64.000   74.000    FALSE 1 1.007   287
# N[5]       58.501  5.098   49.000   58.000   69.000    FALSE 1 1.011   210
# N[6]       61.190  6.112   50.000   61.000   74.000    FALSE 1 1.008   312
# EB[1]      11.218  1.688    7.999   11.163   14.642    FALSE 1 1.000  3000
# EB[2]      13.976  2.161    9.932   13.858   18.382    FALSE 1 1.000  3000
# EB[3]      16.654  2.606   11.743   16.560   21.847    FALSE 1 1.000  3000
# EB[4]      15.085  2.516   10.371   14.963   20.281    FALSE 1 1.000  3000
# EB[5]      13.708  2.339    9.453   13.563   18.695    FALSE 1 1.001  2231
# EB[6]      14.377  2.763    9.430   14.140   20.296    FALSE 1 1.001  2915
# b[1]        0.074  0.012    0.053    0.073    0.098    FALSE 1 1.000  3000
# b[2]        0.104  0.017    0.073    0.104    0.141    FALSE 1 1.000  3000
# b[3]        0.146  0.025    0.101    0.145    0.199    FALSE 1 1.000  3000
# b[4]        0.147  0.029    0.096    0.144    0.212    FALSE 1 1.001  2379
# b[5]        0.145  0.031    0.093    0.143    0.214    FALSE 1 1.002  1177
# b[6]        0.181  0.050    0.103    0.174    0.298    FALSE 1 1.002  1710
# A[1]      151.983  4.074  143.000  152.000  159.000    FALSE 1 1.003   665
# ... [output truncated] ...
# A[6]       81.080  8.215   64.000   81.000   96.000    FALSE 1 1.002  1000
# s[1,1]    624.062  0.791  622.593  624.058  625.670    FALSE 1 1.001  3000
# s[2,1]    629.740  1.956  625.371  629.806  633.406    FALSE 1 1.003   718
# ... [output truncated] ...
# s[199,2] 3401.807 11.931 3381.122 3403.032 3419.268    FALSE 1 1.000  3000
# s[200,2] 3401.362 11.748 3380.839 3402.505 3419.104    FALSE 1 1.000  3000
# zi[1]       1.000  0.000    1.000    1.000    1.000    FALSE 1    NA     1
# ... [output truncated] ...
# zi[200]     0.313  0.464    0.000    0.000    1.000     TRUE 1 1.000  3000
# N.ever    118.920  8.215  104.000  119.000  136.000    FALSE 1 1.002  1000
# z[1,1]      1.000  0.000    1.000    1.000    1.000    FALSE 1    NA     1
# z[2,1]      1.000  0.000    1.000    1.000    1.000    FALSE 1    NA     1
# ... [output truncated] ...
# z[199,6]    0.205  0.404    0.000    0.000    1.000     TRUE 1 1.001  2279
# z[200,6]    0.182  0.386    0.000    0.000    1.000     TRUE 1 1.000  3000


# ~~~~ Fig. 21.3 ~~~~

plot(y=out1$mean$N, x=1:6, ylim=c(0,100), type="b", pch=16, axes=FALSE,
    ylab="Population size (N)", xlab=NA)
segments(1:6, out1$q2.5$N, 1:6, out1$q97.5$N)
axis(1, at=1:6, labels = 2007:2012)
axis(2, las=1)

# ~~~~ Fig. 21.4 ~~~~

library(scales)
co <- viridis_pal(option='E')(20)[c(5, 16)]

plot(0, ylim=c(0, 13), xlim=c(0, 1), frame=FALSE, pch=NA,
    ylab='Posterior density', xlab='Probability', las=1)
lines(density(out1$sims.list$phi), col=co[1], lwd=2)
lines(density(out1$sims.list$gamma), col=co[2] , lwd=2)
legend('topleft', col=co, lwd=rep(2, 2), bty='n',
    legend=c(expression(paste('Adult survival (', phi, ')')),
        expression(paste('Per capita recruitment (', gamma, ')'))))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Make Fig. 21.5
# Define the grid
xgrid <- seq(xlims[1], xlims[2], length.out=61)
ygrid <- seq(ylims[1], ylims[2], length.out=81)

# Grab the MCMC draws for the x and y coordinates of all activity centres (ac) and of z
ac.x <- out1$sims.list$s[,,1]
ac.y <- out1$sims.list$s[,,2]
z <- out1$sims.list$z

library(scales)
library(AHMbook)
cl <- terrain.colors(100, rev=TRUE)[c(1, 1:5*10, 51:100)]
op <- par(mar=c(4, 4, 2, 4), mfrow=c(3, 2))
year <- 2007:2012
# Define density classes
h <- seq(0, 0.3, length.out=100)
for (t in 1:6){
  # Allocate the locations to the defined grid cells
  sx <- cut(ac.x[z[,,t]==1], breaks=xgrid, include.lowest=TRUE)
  sy <- cut(ac.y[z[,,t]==1], breaks=ygrid, include.lowest=TRUE)

  # Cross-tabulation and down-scaling by the number of MCMC draws
  dens <- table(sx, sy) / out1$mcmc.info$n.samples
  # Draw map

  image(x=xgrid, y=ygrid, z=dens, col=cl, asp=1, axes=FALSE, xlim=xlims, ylim=ylims,
      xlab='x coordinate', ylab='y coordinate', main=year[t], zlim=c(0, 0.3))
  axis(1); axis(2, las=1)
}
image_scale(c(0, 0.3), col=cl[seq(1, length(cl), length=13)], digits=2, labels ="breaks",
    cex.legend=0.8)
par(op)

out1$mean$N[3] / 1200 * 100                       # Posterior mean
# > 5.92

out1$q2.5$N[3] / 1200 * 100                       # Lower credible limit
# > 5.42

out1$q97.5$N[3] / 1200 * 100                      # Upper credible limit
# > 6.50


# Make Fig. 21.6
op <- par(mar=c(4, 4, 2, 4), las=1)
# To compute the geometric mean, we only need the density in the first and the last year
# Allocate the locations to the defined grid cells in year 2007
sx <- cut(ac.x[z[,,1]==1], breaks=xgrid, include.lowest=TRUE)
sy <- cut(ac.y[z[,,1]==1], breaks=ygrid, include.lowest=TRUE)
# Cross-tabulation
dens07 <- table(sx, sy) + 0.001

# Allocate the locations to the defined grid cells in year 2012
sx <- cut(ac.x[z[,,6]==1], breaks=xgrid, include.lowest=TRUE)
sy <- cut(ac.y[z[,,6]==1], breaks=ygrid, include.lowest=TRUE)
# Cross-tabulation
dens12 <- table(sx, sy) + 0.001

# Compute geometric mean for each grid
lam <- matrix(NA, nrow=60, ncol=80)
for (i in 1:60){
  for (j in 1:80){
    lam[i,j] <- exp((log(dens12[i,j]) - log(dens07[i,j])) / 5)
  } #j
} #i

# Define the limits of the scale
h <- c(0, 0.9, 0.95, 1, 1.05, 1.10, 20)

# Define colors for that scale
cl <- c(alpha('red', c(0.9, 0.6, 0.3)), alpha('green', c(0.3, 0.6, 0.9)))

# Allocate the actual geom. means into categories
glam <- matrix(.bincode(lam, h), nrow=60, ncol=80)

# Draw map
image(x=xgrid, y=ygrid, z=glam, col=cl, asp=1, axes=FALSE, xlim=xlims, ylim=ylims, xlab='x coordinate',
    ylab='y coordinate', zlim=c(1, 6))
axis(1); axis(2)
par(op)
