# Schaub & Kéry (2022) Integrated Population Models
# Chapter 2 : Bayesian statistical modeling using JAGS
# ----------------------------------------------------

library(IPMbook) ; library(jagsUI)


# 2.7 Using JAGS to fit simple statistical models from R: GLMs and GLMMs
# ======================================================================

# 2.7.7 Generalized linear models with Gaussian random effects (GLMMs)
# --------------------------------------------------------------------

# Set constants in simulation
npop <- 10                                              # Number of populations
grand.mean <- 200                                       # Grand mean of mass
sigma.pop <- 10                                         # SD of population variability in snake mass
beta <- -30                                             # Slope of mass on elevation
sigma.snake <- 50                                       # SD of individual variability of snake mass

# Simulate new snake mass data
set.seed(39)

# Simulate population-level data
mean.mass <- round(rnorm(npop, grand.mean, sigma.pop))  # Average mass (g) in each population at
                                                        # average elevation
elevPop <- sort(runif(npop, 200, 2000))                 # 200 m to 2000 m elevation
elevPopSc <- as.numeric(scale(elevPop))
mean.pop.mass <- mean.mass + beta * elevPopSc
plot(mean.pop.mass, elevPop, pch=16)                    # Not shown: population-level plot

# Simulate snake-level data
nsnakes <- rpois(npop, 20)                              # Number of snakes caught per population
mass <- numeric()
for (j in 1:npop){
  massPop <- rnorm(nsnakes[j], mean.pop.mass[j], sigma.snake)
  mass <- c(mass, massPop)
}

# Create individual-level elevation covariate and population factor
elev <- rep(elevPopSc, nsnakes)
pop <- rep(1:npop, nsnakes)
head(cbind(pop, elev, mass))

# ~~~~ code to plot mass (Fig. 2.18) ~~~~
op <- par(las=1, mar=c(4.1, 4.5, 2, 2), cex=1.5)
plot(mass ~ rep(elevPop, nsnakes), axes=FALSE, pch=16, col=rgb(0,0,0,0.4), xlab='Elevation (m)',
    ylab="Mass (g)", xlim=c(200, 1800), ylim=c(0, 350))
axis(1)
axis(1, at=seq(200, 1700, by=100), tcl=-0.25, labels=NA)
axis(2)
for (i in 1:10){
  segments(elevPop[i]-30, mean.pop.mass[i], elevPop[i]+30, mean.pop.mass[i], col='red', lwd= 2)
}
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Bundle data
jags.data <- list(y=mass, pop=pop, elev=elev, npop=npop, n=length(mass))
str(jags.data)
# List of 5
# $ y   : num [1:178] 214 155 203 214 348 ...
# $ pop : int [1:178] 1 1 1 1 1 1 1 1 1 1 ...
# $ elev: num [1:178] -1.45 -1.45 -1.45 -1.45 -1.45 ...
# $ npop: num 10
# $ n   : int 178

# Write JAGS model file
cat(file="model8.txt", "
model {
  # Priors and linear models
  for (j in 1:npop){
    alpha[j] ~ dnorm(mu.alpha, tau.alpha)               # Random effects
  }
  # Priors for hyper-parameters
  mu.alpha ~ dnorm(0, 1.0E-06)                          # Hyperprior for mean hyperparam
  tau.alpha <- pow(sd.alpha, -2)
  sd.alpha ~ dunif(0, 100)                              # Hyperprior for sd hyperparam

  # Other priors
  beta ~ dnorm(0, 1.0E-06)                              # Slope of mass on elevation
  tau <- pow(sd, -2)
  sd ~ dunif(0, 1000)                                   # 1/residual variance

  # 'Likelihood'
  for (i in 1:n){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha[pop[i]] + beta * elev[i]
  }
}
")

# Initial values
inits <- function() list(mu.alpha=rnorm(1), sd.alpha=runif(1))

# Parameters monitored
parameters <- c("mu.alpha", "sd.alpha", "alpha", "beta", "sd")

# MCMC settings
ni <- 110000; nb <- 10000; nc <- 3; nt <- 2; na <- 1000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out9 <- jags(jags.data, inits, parameters, "model8.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out9) # Not shown
print(out9, 2)

#              mean    sd    2.5%     50%   97.5% overlap0 f Rhat  n.eff
# mu.alpha   195.69  5.64  184.40  195.73  206.87    FALSE 1    1 150000
# sd.alpha    11.66  7.12    0.89   10.80   28.42    FALSE 1    1   6323
# alpha[1]   187.51 10.80  163.20  188.89  205.19    FALSE 1    1  23383
# alpha[2]   203.56  9.62  187.42  202.46  224.92    FALSE 1    1  78358
# alpha[3]   195.27  9.04  176.62  195.43  213.52    FALSE 1    1 145667
# alpha[4]   195.69  8.26  179.01  195.72  212.43    FALSE 1    1 150000
# alpha[5]   201.34  9.49  184.55  200.36  222.35    FALSE 1    1  93368
# alpha[6]   196.33  8.31  179.67  196.26  213.38    FALSE 1    1  90852
# alpha[7]   191.89  8.45  173.73  192.50  207.45    FALSE 1    1 150000
# alpha[8]   183.99 11.19  159.23  185.23  201.78    FALSE 1    1  42527
# alpha[9]   201.49  9.24  185.15  200.59  221.76    FALSE 1    1 130557
# alpha[10]  199.74  9.79  181.62  198.99  221.20    FALSE 1    1 150000
# beta       -30.76  5.83  -42.36  -30.75  -19.18    FALSE 1    1 150000
# sd          47.80  2.62   43.03   47.68   53.28    FALSE 1    1 102355

# Compare the population mean estimates with the known truth
mean.mass
# [1] 198 206 223 182 205 201 201 180 187 207

# Get frequentist restricted MLEs for comparison
library(lme4)
fm <- lmer(mass ~ elev + (1 | as.factor(pop)), REML=TRUE)   # Fit model
summary(fm)                                         # Show summary
ranef(fm)                                           # Give estimates of population mean deviations (not shown)

# Random effects:
#  Groups         Name        Variance Std.Dev.
#  as.factor(Pop) (Intercept)   99.15   9.958
#  Residual                   2242.97  47.360
# Number of obs: 178, groups:  as.factor(Pop), 10

# Fixed effects:
#             Estimate Std. Error t value
# (Intercept)  195.664      4.767  41.048
# elev         -30.777      4.933  -6.239
