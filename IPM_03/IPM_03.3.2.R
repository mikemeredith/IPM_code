# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

library(IPMbook)

  # 3.3 Classical analysis of a matrix population model
# ===================================================

# 3.3.2 Analysis of a matrix population model with parameter uncertainty
# ----------------------------------------------------------------------

# Define mean and SD of the demographic parameters
mean.sj <- 0.3                  # Point estimate of juv. survival
se.sj.e <- 0.03                 # Uncertainty of juv. survival as SE on natural scale
mean.sa <- 0.55                 # Point estimate of ad. survival
se.sa.e <- 0.015                # Uncertainty of ad. survival as SE on natural scale
mean.f1 <- 1.3                  # Point estimate of productivity of 1y females
se.f1.e <- 0.3                  # Uncertainty of 1y productivity as SE on natural scale
mean.fa <- 1.8                  # Point estimate of productivity of ad. females
se.fa.e <- 0.1                  # Uncertainty of ad. productivity as SE on natural scale

# Define number of simulations, vectors and matrices to store results
nsim <- 100000
lambda <- R0 <- GT <- numeric(nsim)
stable.stage <- matrix(NA, ncol=2, nrow=nsim)
sensitivity <- matrix(NA, ncol=4, nrow=nsim)

# Generate demographic values from beta and normal distributions
library(IPMbook)
sj.sim <- rbeta2(nsim, mean.sj, se.sj.e)
sa.sim <- rbeta2(nsim, mean.sa, se.sa.e)
f1.sim <- rnorm(nsim, mean.f1, se.f1.e)
fa.sim <- rnorm(nsim, mean.fa, se.fa.e)

# Perform Monte Carlo simulations
for (s in 1:nsim){
  if(s %% 1000 == 0) {cat(paste("*** Simrep", s, "***\n"))}   # Counter

  # Transition matrix
  A <- matrix(c(sj.sim[s] * f1.sim[s], sj.sim[s] * fa.sim[s], sa.sim[s], sa.sim[s]),
  ncol=2, byrow=TRUE)
  eigenA <- eigen(A)

  # Asymptotic population growth rate
  lambda[s] <- max(Re(eigenA$values))

  # Stable stage distribution
  u <- which.max(Re(eigenA$values))
  revec <- Re(eigenA$vectors[,u])
  stable.stage[s,] <- revec / sum(revec)

  # Lower level sensitivities
  levec <- Re(solve(eigenA$vectors)[u,])
  senmat <- levec %*% t(revec)
  derivmat <- matrix(c(f1.sim[s], fa.sim[s], 0, 0), ncol=2, byrow=TRUE)
  sensitivity[s,1] <- sum(senmat * derivmat)
  derivmat <- matrix(c(sj.sim[s], 0, 0, 0), ncol=2, byrow=TRUE)
  sensitivity[s,2] <- sum(senmat * derivmat)
  derivmat <- matrix(c(0, sj.sim[s], 0, 0), ncol=2, byrow=TRUE)
  sensitivity[s,3] <- sum(senmat * derivmat)
  derivmat <- matrix(c(0, 0, 1, 1), ncol=2, byrow=TRUE)
  sensitivity[s,4] <- sum(senmat * derivmat)

  # Net reproductive rate
  i <- 1:100
  R0[s] <- sj.sim[s] * f1.sim[s] + sj.sim[s] * fa.sim[s] * sum(sa.sim[s]^i)

  # Generation time
  GT[s] <- log(R0[s]) / log(lambda[s])
}

# ~~~ save the results for use in section 3.4.2 ~~~
save(lambda, stable.stage, sensitivity, GT, file="IPM_03.3.2_output.RData")

# ~~~~ code for Table 3.1 ~~~~
# Summaries of the quantities of interest
summary_table <- rbind(
  "Population growth rate"=c(mean(lambda),sd(lambda), quantile(lambda, c(0.025, 0.975))),
      "Stable stage distribution (1y)"=c(mean(stable.stage[,1]), sd(stable.stage[,1]),
      quantile(stable.stage[,1], c(0.025, 0.975))),
  "Stable stage distribution (adults)"=c(mean(stable.stage[,2]), sd(stable.stage[,2]),
      quantile(stable.stage[,2], c(0.025, 0.975))),
  "Growth rate sensitivity of sj"=c(mean(sensitivity[,1]), sd(sensitivity[,1]),
      quantile(sensitivity[,1], c(0.025, 0.975))),
  "Growth rate sensitivity of f1"=c(mean(sensitivity[,2]), sd(sensitivity[,2]),
      quantile(sensitivity[,2], c(0.025, 0.975))),
  "Growth rate sensitivity of fa"=c(mean(sensitivity[,3]), sd(sensitivity[,3]),
      quantile(sensitivity[,3], c(0.025, 0.975))),
  "Growth rate sensitivity of sa"=c(mean(sensitivity[,4]), sd(sensitivity[,4]),
      quantile(sensitivity[,4], c(0.025, 0.975))),
  "Net reproductive rate"=c(mean(R0), sd(R0), quantile(R0, c(0.025, 0.975))),
  "Generation time"=c(mean(GT), sd(GT), quantile(GT, c(0.025, 0.975))))
colnames(summary_table)[1:2] <- c("Mean", "Monte Carlo SE")
round(summary_table, 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code to produce figure 3.11 ~~~~
op <- par(mar=c(4, 4, 3, 0), las=1, cex=1.1, "mfrow")
layout(matrix(1:4, 2, 2, byrow=TRUE), widths=c(1.1, 1), heights=c(1, 1), TRUE)
a <- hist(lambda, nclass=50, col="dodgerblue", main="",
    xlab=expression(paste("Asymptotic population growth rate (", lambda, ")")), prob=TRUE)
mtext("A", at=a$mids[2], cex=1.5)
par(mar=c(4, 2, 3, 2))
a <- hist(stable.stage[,1], nclass=50, col="dodgerblue", main="",
    xlab="Proportion of 1-year old individuals", prob=TRUE)
mtext("B", at=a$mids[2], cex=1.5)
par(mar=c(4, 4, 3, 0))
a <- hist(sensitivity[,1], nclass=50, col="dodgerblue", main="",
    xlab=expression('Sensitivity of '*lambda*' to variation of '*italic(s)[italic(j)]),
    prob=TRUE)
mtext("C", at=a$mids[2], cex=1.5)
par(mar=c(4, 2, 3, 2))
a <- hist(GT, nclass=30, col="dodgerblue", main="", xlab="Generation time", prob=TRUE)
mtext("D", at=a$mids[2], cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
