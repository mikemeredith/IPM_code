# Schaub & Kéry (2022) Integrated Population Models
# Chapter 20 : Kestrel
# --------------------

# Run time for test script 2 hrs, full run 3.5 days!

# 20.4 Component data likelihoods
# ===============================

library(IPMbook); library(jagsUI)
data(kestrel)
str(kestrel)
# List of 4
# $ landData :'data.frame': 15734 obs. of 4 variables:
# ..$ x        : num [1:15734] 920 920 921 921 921 ...
# ..$ y        : num [1:15734] 86.3 87.3 85.3 86.3 87.3 ...
# ..$ elevation: int [1:15734] 1270 1283 1371 1342 1384 1292 1265 1219 ...
# ..$ lakes    : num [1:15734] 0 0 0 0 0 0 0 0 0 0 ...
# $ mhbData  :List of 3
# ..$ x     : num [1:102] 929 935 935 947 947 ...
# ..$ y     : num [1:102] 103.3 95.3 111.3 95.3 111.3 ...
# ..$ tcount: int [1:102, 1:19] NA 0 0 NA NA 0 0 0 0 0 ...
# .. ..- attr(*, "dimnames")=List of 2
# .. .. ..$ : NULL
# .. .. ..$ : chr [1:19] "1999" "2000" "2001" "2002" ...
# $ atlasData:List of 4
# ..$ x    : num [1:574] 923 929 947 953 953 ...
# ..$ y    : num [1:574] 87.3 87.3 119.3 103.3 127.3 ...
# ..$ year : num [1:574] 2015 2013 2013 2013 2015 ...
# ..$ count: num [1:574, 1:3] 0 0 0 1 1 0 0 0 0 2 ...
# .. ..- attr(*, "dimnames")=List of 2
# .. .. ..$ : NULL
# .. .. ..$ : chr [1:3] "count1" "count2" "count3"
# $ drData   :List of 3
# ..$ site        : num [1:24561, 1:3] 1 1 1 1 1 1 1 1 1 1 ...
# .. ..- attr(*, "dimnames")=List of 2
# .. .. ..$ : NULL
# .. .. ..$ : chr [1:3] "site" "x" "y"
# ..$ deadrecovery: num [1:24561, 1:15] 0 0 0 0 0 0 0 0 0 0 ...
# .. ..- attr(*, "dimnames")=List of 2
# .. .. ..$ : NULL
# .. .. ..$ : chr [1:15] "2002" "2003" "2004" "2005" ...
# ..$ age         : chr [1:24561] "juv" "juv" "juv" "juv" ...


# 20.4.1 MHB population count data
# --------------------------------

# Grab MHB total count data as a matrix
C <- as.matrix(kestrel$mhbData$tcount)
dimnames(C) <- NULL

# Produce 'stretched counts', filled up with NAs for unsurveyed quadrats and where the MHB counts
# are filled in at the right place
mhb.coord <- paste(kestrel$mhbData$x, kestrel$mhbData$y, sep=".")
lscape.coord <- paste(kestrel$landData$x, kestrel$landData$y, sep=".")
C.mhb <- array(NA, dim=c(length(kestrel$landData$x), ncol(C)))
mhb.positions <- pmatch(mhb.coord, lscape.coord)
for (i in 1:nrow(C)){
  C.mhb[mhb.positions[i],] <- C[i,]
}


# 20.4.2 Atlas population count data
# ----------------------------------

C.atlas <- kestrel$atlasData$count
dimnames(C.atlas) <- NULL

# Align the atlas data with the modeled domain by determining their position within the landscape
atlas.coord <- paste(kestrel$atlasData$x, kestrel$atlasData$y, sep=".")
atlas.positions <- pmatch(atlas.coord, lscape.coord)


# 20.4.3 Dead-recovery data
# -------------------------

# Check which birds (age) were marked at which site
table(kestrel$drData$site[,1], kestrel$drData$age) # Not shown

# Compute the m-arrays, one for each site and age at ringing
# 1. For individuals ringed as adults (from site 3 only)
a <- with(kestrel$drData, which(age=="ad" & site==3))
ma <- with(kestrel$drData, marrayDead(deadrecovery[a,]))

# 2. For individuals ringed as juveniles (all 7 sites)
mj <- array(0, dim=c(7, 14, 15))
for (i in 1:7){
  a <- with(kestrel$drData, which(age=="juv" & site[,1]==i))
  mj[i,,] <- with(kestrel$drData, marrayDead(deadrecovery[a,]))
}

# Compute the annual number of marked individuals
rel.a <- rowSums(ma)
rel.j <- apply(mj, 2:1, sum)

# Get position of each ringing area within modeled landscape
dr.sites <- unique(kestrel$drData$site)
dr.site.coord <- paste(dr.sites[,"x"], dr.sites[,"y"], sep=".")
dr.positions <- pmatch(dr.site.coord, lscape.coord)
nsites.dr <- length(dr.positions)


# 20.4.4 Basis function approach to the modeling of spatial autocorrelation
# -------------------------------------------------------------------------

library(AHMbook)
library(fields)
coordgrid <- cbind(kestrel$landData$x, kestrel$landData$y)

# Run a knot design
nknots <- 50                                              # This is K in the algebra above
set.seed(50)
knots <- cover.design(R=coordgrid, nd=nknots, nruns=10, num.nn=200,
    max.loop=20)                                          # Takes about 25 min

# Plot the knot design (Fig. 20.4)
library(raster)
r <- rasterFromXYZ(data.frame(x=coordgrid[,1], y=coordgrid[,2], z=kestrel$landData$elevation))
plot(r, col=terrain.colors(20), axes=FALSE, box=FALSE, main="", asp=1)
points(knots$design, col="black", pch=16, cex=1.5)

# Define the Z matrix for the random effects/knot coefficients
knotlocs <- knots$design
omega <- (e2dist(knotlocs, knotlocs)/100)^3               # Scale distances
svd.omega <- svd(omega)
sqrt.omega <- t(svd.omega$v %*% (t(svd.omega$u) * sqrt(svd.omega$d)))
Zk <- (e2dist(coordgrid, knotlocs) / 100)^3               # Scale distances
Zmat <- t(solve(sqrt.omega, t(Zk)))

# Visualize basis vectors
head(Zmat)                                                # Look at first couple of values

# par(mfrow=c(3, 3), mar=c(2, 2, 3, 8), cex.main=1.5, ask=TRUE)
op <- par(mfrow=c(3, 3), mar=c(2, 2, 3, 8), cex.main=1.5,
    ask=dev.interactive(orNone=TRUE))                     # ~~~ only ask if plotting is on-screen
for (i in 1:nknots){
  r <- rasterFromXYZ(data.frame(x=coordgrid[,1], y=coordgrid[,2], z=Zmat[,i]))
  plot(r, col=topo.colors(20), axes=FALSE, box=FALSE, main=paste("Basis vector", i),
      legend=TRUE,axis.args=list(cex.axis=1.5), legend.width=1.5)
  points(knots$design[i,1], knots$design[i,2], col='black', pch=16, cex=1.5)
}
par(op)

# Get values of the Zmat for the atlas data (n = 574 sites)
str(Zmat.at <- Zmat[atlas.positions,])

# Get values of the Zmat for the dead-recovery data (n = 7 sites)
str(Zmat.dr <- Zmat[dr.positions,])


# 20.4.5 Scaling the modelled population size to the nominal 1 km2 area
# ---------------------------------------------------------------------

land.quad <- which(kestrel$landData$lakes < 0.5)


# 20.5 The integrated population model
# ====================================

# Bundle data and produce data overview
jags.data <- list(nsites.region=nrow(C.mhb), nyears.mhb=ncol(C.mhb), C.mhb=C.mhb, nknots=nknots,
    Zmat=Zmat, land.quad=land.quad, nsites.atlas=nrow(C.atlas), nsurveys.atlas=ncol(C.atlas),
    C.atlas=C.atlas, year.atlas=kestrel$atlasData$year, Zmat.at=Zmat.at,nsites.dr=nsites.dr,
    nyears.dr=ncol(kestrel$drData$deadrecovery), mj=mj, rel.j=rel.j, ma=ma,
    rel.a=rel.a, Zmat.dr=Zmat.dr)
str(jags.data)
# List of 18
# $ nsites.region : int 15734
# $ nyears.mhb    : int 19
# $ C.mhb         : int [1:15734, 1:19] NA NA NA NA NA NA NA NA NA NA ...
# $ nknots        : num 50
# $ Zmat          : num [1:15734, 1:50] -0.0666 -0.0647 -0.0654 -0.0635 ...
# $ land.quad     : int [1:15055] 1 2 3 4 5 6 7 8 9 10 ...
# $ nsites.atlas  : int 574
# $ nsurveys.atlas: int 3
# $ C.atlas       : num [1:574, 1:3] 0 0 0 1 1 0 0 0 0 2 ...
# $ year.atlas    : num [1:574] 2015 2013 2013 2013 2015 ...
# $ Zmat.at       : num [1:574, 1:50] -0.05587 -0.04094 0.00386 0.0013 ...
# $ nsites.dr     : int 7
# $ nyears.dr     : int 15
# $ mj            : num [1:7, 1:14, 1:15] 0 0 4 1 4 0 0 0 0 0 ...
# $ rel.j         : num [1:14, 1:7] 22 37 50 152 148 235 117 163 289 318 ...
# $ ma            : num [1:14, 1:15] 0 0 0 0 0 0 0 0 0 0 ...
# $ rel.a         : num [1:14] 69 74 86 64 81 89 64 86 77 45 ...
# $ Zmat.dr       : num [1:7, 1:50] 0.37579 0.42542 0.01412 -0.00927...

# Explore shape and location of some distributions for priors (Fig. 20.6)
op <- par(mfrow=c(2, 2), mar=c(6, 6, 5, 2), cex.axis=1.5, cex.lab=1.5, cex.main=1.5, las=1)
curve(dnorm(x, mean=0, sd=sqrt(1 / 0.1)), 0, 10, col="blue", lwd=3, xlab="", ylab="Density",
    frame=FALSE)
mtext(expression(paste("Prior for intercept of ", lambda)), side=3, line=0)
curve(dbeta(x, shape1=6, shape2=4), 0, 1, col="blue", lwd=3, xlab="", ylab="Density", frame=FALSE)
mtext(expression(paste("Prior for intercept of ", phi)), side=3, line=0)
curve(dnorm(x, mean=0, sd=sqrt(1 / 1)), 0, 3, col="blue", lwd=3, xlab="", ylab="Density", frame=FALSE)
mtext(expression(paste("Prior for intercept of ", gamma)), side=3, line=0)
curve(dbeta(x, shape1=1, shape2=1), 0, 1, col="blue", lwd=3, xlab="", ylab="Density", frame=FALSE)
mtext(expression(paste("Prior for constant ", italic(p))), side=3, line=0)
par(op)


# Write JAGS model file
cat(file = "model1.txt","
model {
  # Priors and linear models
  lambda.int ~ dnorm(0, 0.1)I(0,)                         # Intercept of initial site-specific abundance
  alpha.lam <- log(lambda.int)
  phi.int ~ dbeta(6, 4)                                   # Intercept of survival
  alpha.phi <- logit(phi.int)
  gamma.int ~ dnorm(0, 1)I(0,)                            # Intercept of recruitment
  alpha.gam <- log(gamma.int)
  p ~ dunif(0, 1)                                         # Constant detection prob. (per survey)

  # Fixed parameter for effective survey quadrat size (in km2)
  # with uniform prior for automatic sensitivity analysis/error propagation
  A ~ dunif(0.8, 4)                                       # Effective sampling area of survey

  # Priors on random effects parameters representing the splines
  for (k in 1:nknots){
    b.lam[k] ~ dnorm(0, tau.b.lam)
    b.phi[k] ~ dnorm(0, tau.b.phi)
    b.gam[k] ~ dnorm(0, tau.b.gam)
  }

  # Priors on random effects dispersion (variance among values at knots)
  tau.b.lam <- pow(sd.b.lam, -2)
  tau.b.phi <- pow(sd.b.phi, -2)
  tau.b.gam <- pow(sd.b.gam, -2)
  sd.b.lam ~ dnorm(0, 0.01)I(0,)
  sd.b.phi ~ dnorm(0, 0.01)I(0,)
  sd.b.gam ~ dnorm(0, 0.1)I(0,)

  # Linear models for lambda, phi and gamma for MHB sites
  for (i in 1:nsites.region){                             # Loop over all sites in the modeled domain,
                                                          # including the 102 sites with MHB counts
    lambda[i] <- exp(alpha.lam + smooth.lam[i])
    phi[i] <- ilogit(alpha.phi + smooth.phi[i])
    gamma[i] <- exp(alpha.gam + smooth.gam[i])
    smooth.lam[i] <- inprod(Zmat[i,], b.lam[])            # Random field in lambda
    smooth.phi[i] <- inprod(Zmat[i,], b.phi[])            # Random field in phi
    smooth.gam[i] <- inprod(Zmat[i,], b.gam[])            # Random field in gamma
  }

  # Calculate lambda, phi and gamma for atlas sites
  for (i in 1:nsites.atlas){                              # Loop over 574 sites in the atlas data set
    lamA[i] <- exp(alpha.lam + inprod(Zmat.at[i,], b.lam[]))
    phiA[i] <- ilogit(alpha.phi + inprod(Zmat.at[i,], b.phi[]))
    gamA[i] <- exp(alpha.gam + inprod(Zmat.at[i,], b.gam[]))
  }

  # Linear models for parameters of the dead-recovery model
  for (t in 1:(nyears.dr-1)){
    for (i in 1:nsites.dr){
      logit(sj[i,t]) <- mu.sj + eta.sj[i]
      sa[i,t] <- mean.sa[i]
    } #i
    r[t] <- mean.r
  } #t
  for (i in 1:nsites.dr){
    eta.sj[i] ~ dnorm(0, tau.sj)
    site.sj[i] <- ilogit(mu.sj + eta.sj[i])
    mean.sa[i] <- ilogit(alpha.phi + smooth.lsa[i])       # alpha.phi shared
    smooth.lsa[i] <- inprod(Zmat.dr[i,], b.phi[])         # b.phi also shared
  }

  # Priors for juvenile survival and dead-recovery probability
  mean.sj ~ dunif(0, 1)
  mu.sj <- logit(mean.sj)
  sigma.sj ~ dunif(0, 5)
  tau.sj <- pow(sigma.sj, -2)
  mean.r ~ dunif(0, 1)

  # Longitudinal population count data (MHB, demogr. state-space model)
  for (i in 1:nsites.region){ # Loop over entire region
    # Model for initial condition for MHB data
    N[i,1] ~ dpois(A * lambda[i])
    # Process model over time: our model for population dynamics
    for (t in 1:(nyears.mhb-1)){
      S[i,t+1] ~ dbin(phi[i], N[i,t])                     # Survival process
      R[i,t+1] ~ dpois(gamma[i])                          # 'absolute' recruitment
      N[i,t+1] <- S[i,t+1] + R[i,t+1]
    } #t
  } #i

  # Observation model for the MHB data
  for (i in 1:nsites.region){
    for (t in 1:nyears.mhb){
      C.mhb[i,t] ~ dbin(1 - (1-p)^3, N[i,t])              # Reconcile temporal mismatch
    } #t
  } #i

  # Cross-sectional population count data (atlas, static N-mixture model)
  # Define expected abundance for each year
  for (i in 1:nsites.atlas){                              # Loop over 574 sites in the atlas data set
    # Project expected abundance from 1999 to 2013/2014/2015 (part of calcs)
    for (t in 1:16){
      su[i,t] <- pow(phiA[i], t-1)
    } #t
  } #i

  # Reconcile temporal mismatch between atlas and MHB data
  for (i in 1:182){                                       # Quadrats surveyed in 2013
    lambdaAtlas[i] <- lamA[i] * phiA[i]^14 + gamA[i] * sum(su[i,1:14])
    N.atlas[i] ~ dpois(A * lambdaAtlas[i])
  }
  for (i in 183:387){                                     # Quadrats surveyed in 2014
    lambdaAtlas[i] <- lamA[i] * phiA[i]^15 + gamA[i] * sum(su[i,1:15])
    N.atlas[i] ~ dpois(A * lambdaAtlas[i])
  }
  for (i in 388:574){                                     # Quadrats surveyed in 2015
    lambdaAtlas[i] <- lamA[i] * phiA[i]^16 + gamA[i] * sum(su[i,1:16])
    N.atlas[i] ~ dpois(A * lambdaAtlas[i])
  }

  # Observation model for the atlas data
  for (i in 1:nsites.atlas){                              # Loop over all 574 sites in the atlas data set
    for (j in 1:nsurveys.atlas){                          # Loop over 2-3 occasions
      C.atlas[i,j] ~ dbin(p, N.atlas[i])                  # Assume constant p
    } #j
  } #i

  # Dead-recovery data (dead-recovery model with multinomial likelihood)
  # m-array for adults
  for (t in 1:(nyears.dr-1)){
    ma[t,1:nyears.dr] ~ dmulti(pr.a[t,], rel.a[t])
  }
  for (t in 1:(nyears.dr-1)){
    pr.a[t,t] <- (1-sa[3,t])*r[t]
    for (j in (t+1):(nyears.dr-1)){
      pr.a[t,j] <- prod(sa[3,t:(j-1)]) * (1-sa[3,j]) * r[j]
    } #j
    for (j in 1:(t-1)){
      pr.a[t,j] <- 0
    } #j
  } #t
  for (t in 1:(nyears.dr-1)){
    pr.a[t,nyears.dr] <- 1-sum(pr.a[t,1:(nyears.dr-1)])
  } #t

  # m-array for juveniles (7 sites)
  for (k in 1:nsites.dr){
    for (t in 1:(nyears.dr-1)){
      mj[k,t,1:nyears.dr] ~ dmulti(pr.j[k,t,], rel.j[t,k])
    } #t
    for (t in 1:(nyears.dr-1)){
      pr.j[k,t,t] <- (1-sj[k,t])*r[t]
      for (j in (t+2):(nyears.dr -1)){
        pr.j[k,t,j] <- sj[k,t] * prod(sa[k,(t+1):(j-1)]) * (1-sa[k,j]) * r[j]
      } #j
      for (j in 1:(t-1)){
        pr.j[k,t,j] <- 0
      } #j
    } #t
    for (t in 1:(nyears.dr-2)){
      pr.j[k,t,t+1] <- sj[k,t] * (1-sa[k,t+1]) * r[t+1]
    } #t
    for (t in 1:(nyears.dr-1)){
      pr.j[k,t,nyears.dr] <- 1-sum(pr.j[k,t,1:(nyears.dr-1)])
    } #t
  } #k

  # Derived quantities
  # Estimate population growth rate for each quadrat based on
  # expected abundances rather than on N (to avoid division by zero)
  for (i in 1:nsites.region){
    for (t in 1:(nyears.mhb-1)){
      g[i,t] <- pow(phi[i], t-1)
    } #t
    gr[i] <- pow((lambda[i] * phi[i]^(nyears.mhb-1) + gamma[i] * sum(g[i,]) /
        lambda[i]), 1 / (nyears.mhb-1))
  }

  # Estimate annual abundance in the entire region (excluding lakes)
  Ntot[1] <- sum(lambda[land.quad])                       # Total population size in NW quarter
  for (t in 2:nyears.mhb){
    Ntot[t] <- sum(N[land.quad,t] / A)
  }
}
")

# Initial values
Nst.atlas <- as.numeric(apply(C.atlas, 1, max, na.rm=TRUE) + 1)
Nst.mhb <- C.mhb + 5
Nst.mhb[, 2:19] <- NA                                     # cols 2:19 of N are deterministic, N <- S + R.
R1 <- C.mhb
R1[,1] <- NA
inits <- function(){list(lambda.int=runif(1, 0.4, 0.6), phi.int=runif(1, 0.8, 0.9),
    gamma.int=runif(1, 0.5, 0.9), p=runif(1, 0.6, 0.8), N=Nst.mhb, R=R1 + 1,
    b.lam=runif(jags.data$nknots, -0.5, 0.5), b.phi=runif(jags.data$nknots, -0.5, 0.5),
    b.gam=runif(jags.data$nknots, -0.5, 0.5), mean.r=runif(1, 0, 0.2), N.atlas=Nst.atlas)}

# Parameters monitored
parameters <-c("lambda.int", "alpha.lam", "phi.int", "alpha.phi", "gamma.int", "alpha.gam", "p",
    "sd.b.lam", "sd.b.phi", "sd.b.gam", "mean.sj", "site.sj", "sigma.sj", "mean.sa", "mean.r", "Ntot",
    "b.lam", "b.phi", "b.gam", "lambda", "phi", "gamma", "gr")

# MCMC settings
# ni <- 150000; nb <- 50000; nc <- 3; nt <- 200; na <- 10000
ni <- 1500; nb <- 500; nc <- 3; nt <- 2; na <- 1000  # ~~~ for testing, 2 hrs

# ~~~ This run requires approx. 5GB of memory. If insufficient memory is available to run all
# the chains, the error message will read 
# "Error in unserialize(node$con) : error reading from connection"

# Call JAGS from R (ART 88 h!) and check convergence
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
# ~~~ Good idea to save it after waiting so long! ~~~
save(out1, file="IPM_20_out1.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
traceplot(out1)

# ~~~~~ Fig. 20.7 ~~~~~
traceplot(out1, c("lambda.int", "alpha.lam", "phi.int", "alpha.phi", "gamma.int",
    "alpha.gam", "sd.b.lam", "sd.b.phi", "sd.b.gam"))
# ~~~~~~~~~~~~~~~~~~~~~


# 20.6 Results
# ============

print(out1$summary[1:46,-c(4,6)], 3)

#                 mean         sd      2.5%       50%      97.5% Rhat n.eff overlap0     f
# lambda.int    0.3771    0.37939    0.1384    0.2763     1.4014 1.01   664        0 1.000
# alpha.lam    -1.1866    0.56794   -1.9779   -1.2862     0.3375 1.01   414        1 0.949
# phi.int       0.7015    0.14452    0.4001    0.7124     0.9297 1.01   172        0 1.000
# alpha.phi     0.9692    0.77920   -0.4051    0.9072     2.5815 1.02   162        1 0.900
# gamma.int     0.0856    0.04414    0.0278    0.0766     0.1967 1.06    82        0 1.000
# alpha.gam    -2.5742    0.49159   -3.5831   -2.5697    -1.6260 1.05    75        0 1.000
# p             0.2777    0.01882    0.2427    0.2766     0.3178 1.01   185        0 1.000
# sd.b.lam      0.5732    0.66389    0.0143    0.3371     2.3572 1.01   450        0 1.000
# sd.b.phi      2.7273    1.05800    0.8968    2.6225     5.1486 1.05    48        0 1.000
# sd.b.gam      0.8178    0.36308    0.3109    0.7451     1.7188 1.03   127        0 1.000
# mean.sj       0.6231    0.08182    0.4486    0.6261     0.7795 1.00  1500        0 1.000
# site.sj[1]    0.6702    0.06186    0.5354    0.6725     0.7809 1.00   674        0 1.000
# site.sj[2]    0.5112    0.04943    0.4132    0.5114     0.6081 1.00  1389        0 1.000
# site.sj[3]    0.8078    0.03296    0.7341    0.8104     0.8615 1.01   388        0 1.000
# site.sj[4]    0.5677    0.07773    0.4115    0.5697     0.7155 1.00  1500        0 1.000
# site.sj[5]    0.5694    0.04736    0.4733    0.5698     0.6597 1.00  1280        0 1.000
# site.sj[6]    0.6458    0.16850    0.2660    0.6616     0.9439 1.00   441        0 1.000
# site.sj[7]    0.5693    0.08679    0.3872    0.5756     0.7224 1.00   680        0 1.000
# sigma.sj      0.7901    0.42492    0.3249    0.6845     1.8401 1.01   602        0 1.000
# mean.sa[1]    0.8952    0.03187    0.8203    0.9000     0.9424 1.00   510        0 1.000
# mean.sa[2]    0.8146    0.04639    0.7120    0.8195     0.8913 1.01   235        0 1.000
# mean.sa[3]    0.9474    0.01518    0.9137    0.9495     0.9666 1.05   207        0 1.000
# mean.sa[4]    0.8309    0.05409    0.7049    0.8374     0.9168 1.00   832        0 1.000
# mean.sa[5]    0.7689    0.04217    0.6840    0.7688     0.8487 1.00   514        0 1.000
# mean.sa[6]    0.9167    0.03380    0.8361    0.9228     0.9648 1.00   363        0 1.000
# mean.sa[7]    0.8736    0.05797    0.7252    0.8861     0.9504 1.05    67        0 1.000
# mean.r        0.0396    0.00383    0.0326    0.0394     0.0477 1.01   289        0 1.000
# Ntot[1]    3757.6250  879.61423 2288.8791 3667.5283  5674.0263 1.00  1500        0 1.000
# Ntot[2]    4586.0053  856.73316 3195.7486 4496.2487  6475.2987 1.00  1500        0 1.000
# Ntot[3]    5236.5288  892.89131 3750.7405 5160.4451  7234.0252 1.00  1500        0 1.000
# Ntot[4]    5773.9687  949.04202 4164.5854 5683.4273  7935.3141 1.00  1500        0 1.000
# Ntot[5]    6229.0835 1014.34675 4483.0932 6124.4378  8545.4490 1.00  1500        0 1.000
# Ntot[6]    6628.2767 1083.24952 4752.5860 6551.7435  9174.1480 1.00  1500        0 1.000
# Ntot[7]    6971.9769 1148.07832 4954.9725 6877.1834  9607.1014 1.00  1500        0 1.000
# Ntot[8]    7274.4552 1211.63479 5162.8650 7165.3406  9999.6642 1.00  1500        0 1.000
# Ntot[9]    7543.5614 1270.96517 5295.5986 7432.9123 10429.5929 1.00  1500        0 1.000
# Ntot[10]   7783.4523 1326.65793 5443.9203 7664.0278 10809.4471 1.00  1500        0 1.000
# Ntot[11]   7997.1026 1378.66310 5574.9320 7883.5742 11130.7045 1.00  1500        0 1.000
# Ntot[12]   8192.6759 1427.22263 5685.0568 8086.6733 11396.9888 1.00  1500        0 1.000
# Ntot[13]   8367.8689 1474.01509 5786.9385 8245.9030 11610.4482 1.00  1500        0 1.000
# Ntot[14]   8532.3876 1517.62460 5887.2966 8416.9704 11879.4728 1.00  1500        0 1.000
# Ntot[15]   8675.6203 1559.33544 5959.8361 8558.4323 12121.7904 1.00  1500        0 1.000
# Ntot[16]   8804.8099 1596.13958 6045.7682 8678.9731 12362.7632 1.00  1500        0 1.000
# Ntot[17]   8923.5676 1630.29313 6109.6724 8808.8073 12619.2903 1.00  1500        0 1.000
# Ntot[18]   9035.4767 1664.85550 6191.8468 8912.2445 12822.9241 1.00  1500        0 1.000
# Ntot[19]   9135.8852 1696.37655 6220.5593 8989.3095 13014.6297 1.00  1500        0 1.000


# Produce Fig. 20.8
op <- par(mfrow=c(2, 3), mar=c(4, 4, 6, 9), cex.main=1.6)

# Random field of initial abundance
r <- rasterFromXYZ(data.frame(x=coordgrid[,1], y=coordgrid[,2], z=out1$mean$lambda))
plot(r, col=topo.colors(20), axes=FALSE, box=FALSE, asp=1, main="Initial density",
    axis.args=list(cex.axis=2), legend.width=3)

# Random field of apparent survival
r <- rasterFromXYZ(data.frame(x=coordgrid[,1], y=coordgrid[,2], z=out1$mean$phi))
plot(r, col=topo.colors(20), axes=FALSE, box=FALSE, asp=1, main="Survival",
    axis.args=list(cex.axis=2), legend.width=3)

# Random field of (absolute) recruitment
r <- rasterFromXYZ(data.frame(x=coordgrid[,1], y=coordgrid[,2], z=out1$mean$gamma))
plot(r, col=topo.colors(20), axes=FALSE, box=FALSE, asp=1, main="Recruitment",
    axis.args=list(cex.axis=2), legend.width=3)

# Uncertainty in the random field of initial abundance (% CV)
r <- rasterFromXYZ(data.frame(x=coordgrid[,1], y=coordgrid[,2], z=100 *
    out1$sd$lambda / out1$mean$lambda))
plot(r, col=topo.colors(20), axes=FALSE, box=FALSE, asp=1, main="Uncertainty of density (% CV)",
    axis.args=list(cex.axis=2), legend.width=3)

# Uncertainty in the random field of survival (% CV)
r <- rasterFromXYZ(data.frame(x=coordgrid[,1], y=coordgrid[,2], z=100 *
    out1$sd$phi /out1$mean$phi))
plot(r, col=topo.colors(20), axes=FALSE, box=FALSE, asp=1, main="Uncertainty of survival (% CV)",
    axis.args=list(cex.axis=2), legend.width=3)

# Uncertainty in the random field of recruitment (% CV)
r <- rasterFromXYZ(data.frame(x=coordgrid[,1], y=coordgrid[,2], z=100 *
    out1$sd$gamma / out1$mean$gamma))
plot(r, col=topo.colors(20), axes=FALSE, box=FALSE, asp=1, main="Uncertainty of recruitment (% CV)",
    axis.args=list(cex.axis=2), legend.width=3)
par(op)


# ~~~~ Produce Fig. 20.9 ~~~~
op <- par(mfrow=c(1, 2), mar=c(4, 4, 6, 9), cex.main=1.6)
# Random field of population growth rate
r <- rasterFromXYZ(data.frame(x=coordgrid[,1], y=coordgrid[,2], z=out1$mean$gr))
plot(r, col=topo.colors(20), axes=FALSE, box=FALSE, main="Annual population growth rate",
    asp=1, axis.args=list(cex.axis=2), legend.width=2)

# Uncertainty in the random field of population growth rate (% CV)
r <- rasterFromXYZ(data.frame(x=coordgrid[,1], y=coordgrid[,2],
    z=100 * (out1$sd$gr / out1$mean$gr)))
plot(r, col=topo.colors(20), axes=FALSE, box=FALSE, main="Uncertainty of growth rate (% CV)",
    asp=1, axis.args=list(cex.axis=2), legend.width=2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Produce Fig. 20.10 ~~~~
library(scales)
op <- par(mfrow=c(1, 2), mar=c(6, 8, 4, 2), cex.lab=1.5, cex.axis=1.5, las=1)
# Plot of estimated trajectory of the total population size
plot(1999:2017, out1$mean$Ntot, xlab="Year", ylab=NA, pch=16, cex=1.5, axes=FALSE, ylim=c(0, 13000))
mtext("Population size", side=2, line=5, las=0, cex=1.5)
polygon(c(1999:2017, 2017:1999), c(out1$q2.5$Ntot, rev(out1$q97.5$Ntot)),
    border=NA, col=alpha("grey", 0.5))
points(1999:2017, out1$mean$Ntot, pch = 16, type="p", cex=1.5)
axis(1, at=1999:2017, tcl=-0.25, labels=NA)
axis(1, at=c(2000, 2005, 2010, 2015), tcl=-0.5, labels=c(2000, 2005, 2010, 2015))
axis(2)

# Estimate the increase in the population between 1999 and 2017
increase <- out1$sims.list$Ntot[,19] / out1$sims.list$Ntot[,1]
hist(increase, breaks=30, col=alpha("grey", 0.5), main="", xlab="Increase 1999-2017")
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
