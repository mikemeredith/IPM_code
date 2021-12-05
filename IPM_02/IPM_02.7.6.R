# Schaub & Kéry (2022) Integrated Population Models
# Chapter 2 : Bayesian statistical modeling using JAGS
# ----------------------------------------------------

library(IPMbook) ; library(jagsUI)


# 2.7 Using JAGS to fit simple statistical models from R: GLMs and GLMMs
# ======================================================================

# 2.7.6 Normal linear regression or Gaussian generalized linear models
# --------------------------------------------------------------------

# Choose constants and create population indicator
set.seed(88)
npop <- 3                                           # Number of populations
nsnakes <- rpois(npop, 20)                          # Number of snakes caught per population
sigma <- 50                                         # Snake-to-snake variability in mass
Pop <- as.factor(rep(c("A", "B", "C"), nsnakes))    # Create population indicator (factor)

# Generate mass data
mean.mass <- round(runif(npop, 100, 300))           # Pop. average mass (g)
massA <- rnorm(nsnakes[1], mean.mass[1], sigma)     # Mass data in population A
massB <- rnorm(nsnakes[2], mean.mass[2], sigma)     # Mass data in population B
massC <- rnorm(nsnakes[3], mean.mass[3], sigma)     # Mass data in population C
mass <- round(c(massA, massB, massC))               # Round and stack mass data

table(Pop)                                          # Look at sample size per pop. (not shown)
tapply(mass, Pop, mean)                             # Look at sample means per population
mean.mass                                           # Compare with known true means (not shown)
#        A        B        C
# 131.6111 241.5000 246.8333


# ~~~~ extra code for plotting ~~~~
# Plot mass
plot(mass ~ as.numeric(Pop), axes=FALSE, pch=16, cex=2, col=rgb(0, 0, 0, 0.4),
    xlab='', ylab="Mass (g)", xlim=c(0.8, 3.2), ylim=c(0, 350))
axis(1, at=1:3, c("Pop A", "Pop B", "Pop C"), las=1)
axis(2)
segments(0.8, mean.mass[1], 1.2, mean.mass[1], col='red', lwd=2)
segments(1.8, mean.mass[2], 2.2, mean.mass[2], col='red', lwd=2)
segments(2.8, mean.mass[3], 3.2, mean.mass[3], col='red', lwd=2)

# Figure 2.17
s <- as.numeric(table(Pop))
op <- par(las=1, mar=c(4.1, 4.5, 2, 2), cex=1.5)
plot(mass[Pop=='A'], rep(3, s[1]), axes=FALSE, pch=16, ylab="", col=rgb(0,0,0,0.4),
    xlab="Mass (g)", ylim=c(0.8, 3.2), xlim=c(0, 350))
points(mass[Pop=='B'], rep(2, s[2]), pch=16, col=rgb(0,0,0,0.4))
points(mass[Pop=='C'], rep(1, s[3]), pch=16, col=rgb(0,0,0,0.4))
axis(2, at = 1:3, c("Pop C", "Pop B", "Pop A"), las = 1)
axis(1)
segments(mean.mass[3], 0.8, mean.mass[3], 1.2, col='red', lwd=2)
segments(mean.mass[2], 1.8, mean.mass[2], 2.2, col='red', lwd=2)
segments(mean.mass[1], 2.8, mean.mass[1], 3.2, col='red', lwd=2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Bundle data
pop <- as.numeric(Pop)                              # Code factor levels as 1:3
jags.data <- list(y=mass, pop=pop, n=length(mass))
str(jags.data)
# List of 3
# $ mass: num [1:70] 206 139 60 84 118 142 138 175 133 146 ...
# $ pop : num [1:70] 1 1 1 1 1 1 1 1 1 1 ...
# $ n   : int 70

# Write JAGS model file
cat(file="model7.txt", "
model {
  # Priors and linear models
  for (k in 1:3){
    alpha[k] ~ dnorm(0, 1.0E-06)                      # Population means
    # tau[k] <- pow(sigma[k], -2)                     # For heteroscedasticity
    # sigma[k] ~ dunif(0, 1000)
  }
  # Homogeneous residual variance (homoscedasticity)
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 1000)

  # Likelihood for normal linear model
  for (i in 1:n){
    y[i] ~ dnorm(mu[i], tau)                          # Homoscedasticity
    # y[i] ~ dnorm(mu[i], tau[pop[i]])                # Heteroscedasticity
    mu[i] <- alpha[pop[i]]
  }

  # Derived quantities: contrasts between population means
  diff.12 <- alpha[2]-alpha[1]
  diff.13 <- alpha[3]-alpha[1]
  diff.23 <- alpha[3]-alpha[2]
}
")

# Initial values
inits <- function(){list(alpha=rnorm(3, 100, 300))}

# Parameters monitored
parameters <- c("alpha", "sigma", "diff.12", "diff.13", "diff.23")

# MCMC settings
ni <- 60000; nb <- 10000; nc <- 3; nt <- 5; na <- 1000

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out8 <- jags(jags.data, inits,parameters, "model7.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out8, layout=c(2, 2)) # Not shown
print(out8, 2)

#            mean    sd   2.5%    50%  97.5% overlap0    f Rhat n.eff
# alpha[1] 131.60 10.53 111.10 131.48 152.59    FALSE 1.00    1 11010
# alpha[2] 241.58  9.46 222.92 241.56 260.20    FALSE 1.00    1 30000
# alpha[3] 246.80  8.16 230.86 246.79 262.89    FALSE 1.00    1 12920
# sigma     44.60  3.92  37.69  44.32  53.01    FALSE 1.00    1 30000
# diff.12  109.98 14.22  81.80 110.09 137.43    FALSE 1.00    1 13747
# diff.13  115.20 13.27  88.95 115.28 141.07    FALSE 1.00    1  7177
# diff.23    5.22 12.54 -19.51   5.32  30.04     TRUE 0.66    1 30000

# Compare the population mean estimates with the known truth
mean.mass                                           # For population 1, 2 and 3 (or A, B and C)
# [1] 106 252 235

# Get frequentist MLEs for comparison
# summary(lm(mass ~ Pop))                           # Treatment contrast parameterization
summary(lm(mass ~ Pop-1))                           # Means parameterization of factor levels

# Coefficients:
#      Estimate Std. Error t value Pr(>|t|)
# popA  131.611     10.315   12.76   <2e-16 ***
# popB  241.500      9.331   25.88   <2e-16 ***
# popC  246.833      7.990   30.89   <2e-16 ***

# Residual standard error: 43.76 on 67 degrees of freedom
# Multiple R-squared:  0.9639,    Adjusted R-squared:  0.9622
# F-statistic: 595.7 on 3 and 67 DF,  p-value: < 2.2e-16
