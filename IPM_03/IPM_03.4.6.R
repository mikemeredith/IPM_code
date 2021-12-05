# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

library(IPMbook) ; library(jagsUI)

# 3.4 Analysis of matrix population models with Markov Chain
#     Monte Carlo (MCMC) software
# ==========================================================

# 3.4.6 Matrix population models with density-dependence and demographic stochasticity
# ------------------------------------------------------------------------------------

# Define means of the demographic rates and strength of density dependence
mean.sj <- 0.3
mean.sa <- 0.55
f1.int <- 2.3           # Productivity of 1y females when population size is 0
f1.beta <- -0.02        # Strength of density dependence on 1y productivity
fa.int <- 2.3           # Productivity of adult females when population size is 0
fa.beta <- -0.01        # Strength of density dependence on ad productivity

# Define the number of years with predictions
T <- 200

# Bundle data
jags.data <- list(sj=mean.sj, sa=mean.sa, f1.int=f1.int, f1.beta=f1.beta, fa.int=fa.int,
    fa.beta=fa.beta, T=T)

# Write JAGS model file
cat(file="model8.txt", "
model {
  # Model for initial state
  N[1,1] <- 10
  N[2,1] <- 10

  # Loop over time
  for (t in 1:T){
    # Calculate actual productivity
    f1[t] <- f1.int + f1.beta * (N[1,t] + N[2,t])
    fa[t] <- fa.int + fa.beta * (N[1,t] + N[2,t])

    # Population model
    N[1,t+1] ~ dpois(sj * (f1[t] * N[1,t] + fa[t] * N[2,t]))
    N[2,t+1] ~ dbin(sa, (N[1,t] + N[2,t]))
    extinct[t] <- equals(N[1,t+1] + N[2,t+1], 0)          # Determines whether
        # population is still thriving (extinct = 0) or went extinct (extinct = 1)
  }
}
")

# Initial values
N <- matrix(NA, nrow=2, ncol=T+1)
N[,2:(T+1)] <- 10
inits <- function(){list(N=N)}

# Parameters monitored
parameters <- c("N", "f1", "fa", "extinct")

# MCMC settings
ni <- 1000; nt <- 1; nb <- 0; nc <- 1; na <- 0

# Call JAGS (ART <1 min) and summarize results
out8 <- jags(jags.data, inits, parameters, "model8.txt", n.adapt=na, n.chains=nc, n.thin=nt,
    n.iter=ni, n.burnin=nb, DIC=FALSE)
print(out8, 4)

#                 mean     sd    2.5%    50%   97.5% overlap0 f
# N[1,1]       10.0000 0.0000 10.0000 10.000 10.0000    FALSE 1
# N[2,1]       10.0000 0.0000 10.0000 10.000 10.0000    FALSE 1
# N[1,2]       11.9210 3.3838  6.0000 12.000 19.0000    FALSE 1
# N[2,2]       11.0280 2.2804  7.0000 11.000 15.0000    FALSE 1
# N[1,3]       13.4220 4.2102  6.0000 13.000 22.0000    FALSE 1
# N[2,3]       12.6680 3.3665  6.0000 12.000 19.0250    FALSE 1
# [ ... output truncated ... ]
# f1[1]         1.9000 0.0000  1.9000  1.900  1.9000    FALSE 1
# f1[2]         1.8410 0.0814  1.6795  1.840  2.0000    FALSE 1
# f1[3]         1.7782 0.1257  1.5195  1.780  2.0200    FALSE 1
# f1[4]         1.7197 0.1594  1.4000  1.720  2.0200    FALSE 1
# [ ... output truncated ... ]
# fa[1]         2.1000 0.0000  2.1000  2.100  2.1000    FALSE 1
# fa[2]         2.0705 0.0407  1.9897  2.070  2.1500    FALSE 1
# fa[3]         2.0391 0.0628  1.9097  2.040  2.1600    FALSE 1
# fa[4]         2.0099 0.0797  1.8500  2.010  2.1600    FALSE 1
# [ ... output truncated ... ]


# Carrying capacity
mean(rowSums(out8$sims.list$N[,,T+1]))
# [1] 53.382

# Productivity at carrying capacity
mean(out8$sims.list$f1[,200])
# [1] 1.22968
mean(out8$sims.list$fa[,200])
# [1] 1.76484
