# This is the catch all script
# When I'm testing stuff and don't know where it goes yet

# F. Javiera Rudolph
# f.javierarudolph at gmail dot com
# or submit a Github issue


# Function arguments:
# nunits: Number of main units (large quadrats)
# nsubunits: Number of subunits (nested subsamples within
#    each main unit
# nreps: Number of rep surveys in every subunit
# mean.psi: Mean large-scale, unit-level occupancy probability (psi)
# beta.Xpsi: effect on psi of covariate A (at main unit level)

# not using set to 0
# sd.logit.psi: SD of logit(psi), unstructured site variation in psi

# mean.theta: Mean small-scale (subunit) occupancy probability (theta)
# theta.time.range: range of theta 'intercepts' for subunits
# beta.Xtheta: effect on theta of covariate B (at subunit level)

# not using set to 0
# sd.logit.theta: SD of logit(theta), unstructured site variation in theta


# mean.p: Mean per-survey detection probability
# p.time.range: range of p 'intercepts' for replicates
# beta.Xp: effect on p of covariate C (unit by subunit by replicate)

# not using set to 0
# sd.logit.p: SD of logit(p)

nunits = 20; nsubunits = 4; nreps = 3
nstates = 3

# grid cell level
mean.psi = 0.8; beta.Xpsi = 0; sd.logit.psi = 0
mean.r = 0.3; beta.Xr = 0; sd.logit.r = 0

# site level
mean.theta = 0.6; theta.time.range = c(0,0); beta.Xtheta = 0; sd.logit.theta = 0;

theta1 <- 0.65
theta2 <- 0.8
theta3 <- 0.7

mean.theta <- matrix(c(1,0,0,
                  1-theta1, theta1, 0,
                  1-theta2, theta2*(1-theta3), theta2*theta3), nrow = 3, byrow = TRUE)

# Observation
mean.p = 0.4; p.time.range = c(0,0); beta.Xp = 0; sd.logit.p = 0;

p2 <- 0.9
p3 <- c(0.1, 0.2, 0.7)

mean.p<- matrix(c(1,0,0,
                1-p2, p2, 0,
                p3[1], p3[2], p3[3]), nrow = 3, byrow = TRUE)

show.plot = TRUE; verbose = TRUE


# Create data structures
z <- array(NA, dim = nunits)  # Unit occurrence
omega <- array(NA, dim = c(nunits, nstates))
a <- array(NA, dim = c(nunits, nsubunits)) # Subunit
theta <- array(NA, dim = c(nunits, nsubunits, nstates))
y <- array(NA, dim=c(nunits, nsubunits, nreps) ) # Rep


# Create standardised covariate values
covA <- as.numeric(array(runif(nunits, -2, 2), dim = nunits))
covB <- array(runif(nunits*nsubunits, -2, 2),
              dim = c(nunits, nsubunits))
covC <- array(runif(nunits*nsubunits*nreps, -2, 2),
              dim=c(nunits, nsubunits, nreps) )

# Simulate psi, theta and p and plot all
psi <- plogis(qlogis(mean.psi) + beta.Xpsi * covA + rnorm(nunits, 0, sd.logit.psi))
r <- plogis(qlogis(mean.r) + beta.Xr * covA + rnorm(nunits, 0, sd.logit.r))
omega[,1] <- 1-psi
omega[,2] <- psi*(1-r)
omega[,3] <- psi * r

theta.time.effect <- runif(nsubunits, theta.time.range[1], theta.time.range[2])
p.time.effect <- runif(nreps, p.time.range[1], p.time.range[2])

for(i in 1:nunits){
  for(j in 1:nsubunits){
    theta[i,j,] <- plogis(qlogis(mean.theta) + theta.time.effect[j] + (beta.Xtheta*covB)[,j] + array(rnorm(nunits*nsubunits, 0, sd.logit.theta), dim = c(nunits, nsubunits))[,j])
    for(k in 1:nreps){
      p[,j,k] <- plogis(qlogis(mean.p[]) + p.time.effect[k] + (beta.Xp*covC)[,j,k]+ array(rnorm(nunits*nsubunits*nreps, 0,sd.logit.p),dim =c(nunits, nsubunits, nreps))[,j,k])
    }
  }
}


theta <- mean.theta
p <- mean.p


# Sample three nested Bernoulli distributions
# with probabilities psi, z*theta and a * p
for (i in 1:nunits) {
  draw1 <- rmultinom(n = 1, size = 1, prob = omega[i,])
  z[i] <- which(draw1 == 1)
}

for(i in 1:nunits){
  for (j in 1:nsubunits) {
    draw1 <- rmultinom(n = 1, size = 1, prob = theta[z[i],])
    a[i, j] <- which(draw1 == 1)
  }
}

for(i in 1:nunits){
  for(j in 1:nsubunits){
    for (k in 1:nreps) {
      draw1 <- rmultinom(n=1, size = 1, prob = p[a[i,j],])
      y[i,j,k] <- which(draw1 == 1)
    } # survey
  } # subunit
} # unit

sum.z <- sum(z)
sum.z.a <- sum(apply(a, 1, sum)>0)
obs.sum.z <- sum(apply(apply(y, c(1,2), max), 1, max))
if(verbose) {
  cat(" Occupied units:                           ", sum.z, "\n",
      "Units with >=1 occupied, surveyed subunit:", sum.z.a, "\n",
      "Observed number of occupied units:        ", obs.sum.z, "\n",
      "\n")
}

# Visualisation of covariate relationships of psi, theta and p
if(show.plot) {
  op <- par(mfrow = c(1,3), mar = c(5,5,5,2), cex.lab = 1.5, cex.axis = 1.5) ; on.exit(par(op))
  tryPlot <- try( {
    plot(covA, psi, xlab = "Unit covariate A", ylab = "psi", ylim = c(0,1), main = "Large-scale occupancy probability (psi)", frame = FALSE)
    curve(plogis(qlogis(mean.psi) + beta.Xpsi * x), -2, 2, col = "red", lwd = 3, add = TRUE)
    plot(covB, theta, xlab = "Unit-subunit covariate B", ylab = "theta", ylim = c(0,1), main = "Small-scale occupancy probability/availability \n(theta) (red - time variation)", frame = FALSE)
    for(j in 1:nsubunits){
      curve(plogis(qlogis(mean.theta) + theta.time.effect[j] +
                     beta.Xtheta * x), -2, 2, lwd = 2, col = "red", add = T)
    }
    plot(covC, p, xlab = "Unit-subunit-rep covariate C", ylab = "p", ylim = c(0,1), main = "Detection probability (p) \n (red - replicate variation)", frame = FALSE)
    for(k in 1:nreps){
      curve(plogis(qlogis(mean.p) + p.time.effect[k] +
                     beta.Xp * x), -2, 2, lwd = 2, col = "red", add = T)
    }
  }, silent = TRUE)
  if(inherits(tryPlot, "try-error"))
    tryPlotError(tryPlot)
}


# Bundle and summarize data set
str(win.data <- list(y = y, n.cell = dim(y)[1], n.sites = dim(y)[2],
                      n.surveys = dim(y)[3], covA = covA, covB = covB, covC = covC))

# Define model
sink("jags_txt/model.txt")
cat("
model {
  # Priors and model for params
  int.psi ~ dunif(0,1)         # Intercept of occupancy probability
  int.r ~ dunif(0,1)
  for(t in 1:n.sites){
    int.theta[t,1] ~ dunif(0,1) # Intercepts availability probability
    int.theta[t,2] ~ dunif(0,1)
    int.theta[t,3] ~ dunif(0,1)
  }
  for(t in 1:n.surveys){
    int.p[t,1] ~ dunif(0,1)     # Intercepts detection probability
    int.p[t,2] ~ dunif(0,1)
    int.p[t,3] ~ dunif(0,1)
  }
  beta.lpsi ~ dnorm(0, 0.1)    # Slopes of three covariates
  beta.ltheta ~ dnorm(0, 0.1)
  beta.lp ~ dnorm(0, 0.1)
  # 'Likelihood' (or basic model structure)
  for (i in 1:n.cell){
    # Occurrence in pond i
    z[i] ~ dcat(omega[i,])
    logit(psi[i]) <- logit(int.psi) + beta.lpsi * covA[i]
    logit(r[i]) <- logit(int.r) + beta.lr * covA[i]
    omega[i,1] <- 1-psi[i]
    omega[i,2] <- psi[i] * (1-r[i])
    omega[i,3] <- psi[i] * r[i]
    for (j in 1:n.sites){
      # Occurrence in sample j
      a[i,j] ~ dbern(mu.a[i,j])
      mu.a[i,j] <- z[i] * theta[i,j]
      logit(theta[i,j]) <- logit(int.theta[j]) + beta.ltheta * covB[i,j]
      for (k in 1:n.surveys){
        # PCR detection error process in sample k
        y[i,j,k] ~ dbern(mu.y[i,j,k])
        mu.y[i,j,k] <- a[i,j] * p[i,j,k]
        logit(p[i,j,k]) <- logit(int.p[k]) + beta.lp * covC[i,j,k]
      }
    }
    tmp[i] <- step(sum(a[i,])-0.1)
  }
  # Derived quantities
  sum.z <- sum(z[])   # Total # of occupied ponds in sample
  sum.a <- sum(tmp[]) # Total # of ponds with presence in <=1 of the 5 samples
} # end model
",fill=TRUE)
sink()

# Initial values
zst <- apply(y, 1, max)        # inits for presence (z)
ast <- apply(y, c(1,2), max)   # inits for availability (a)
inits <- function() list(z = zst, a = ast, int.psi = 0.5, beta.lpsi = 0)

# Parameters monitored
params <- c("int.psi", "int.r", "int.theta", "int.p", "beta.lpsi", "beta.lr", "beta.ltheta",
            "beta.lp", "sum.z", "sum.a")

# MCMC setting
ni <- 5000   ;   nt <- 2   ;   nb <- 1000   ;   nc <- 3

# Call WinBUGS and summarize posterior
out <- jagsUI::jags(data = win.data, inits = inits, params, "model.txt",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
print(out, 3)
sum.z




