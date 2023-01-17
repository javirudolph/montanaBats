# TEST SIMULATION AND ESTIMATION
# Attempting a dynamic multistate multiscale occupancy model

library(jagsUI)


ncells = 10; nsites = 4; nreps = 5; nyears = 3; nstates = 3

# cell occupancy parameters
psi = 0.4; r = 0.3

omega <- array(NA, dim = c(ncells, nstates))
omega[,1] <- 1-psi
omega[,2] <- psi*(1-r)
omega[,3] <- psi * r

# site occupancy parameters
theta1 = 0.65; theta2 = 0.8; theta3 = 0.7
theta <- matrix(c(1,0,0,
                  1-theta1, theta1, 0,
                  1-theta2, theta2*(1-theta3), theta2*theta3), nrow = 3, byrow = TRUE)

# detection parameters
p2 = 0.9; p3 = c(0.1, 0.2, 0.7)
mean.p<- matrix(c(1,0,0,
                  1-p2, p2, 0,
                  p3[1], p3[2], p3[3]), nrow = 3, byrow = TRUE)

# transition matrix parameters
gamma1 = 0.001; gamma2 = 0.001; phi1 = 0.5; G = 0.01
phi2 = 0.3; D = 0.5

Phi <- matrix(c(1-gamma1, gamma1*(1-gamma2), gamma1*gamma2,
                1-phi1, phi1*(1-G), phi1*G,
                1-phi2, phi2*D, phi2*(1-D)), nrow = 3, byrow = TRUE)


# create empty arrays
z <- array(NA, dim = c(ncells, nyears))
u <- array(NA, dim = c(ncells, nsites, nyears))
y <- array(NA, dim = c(ncells, nsites, nreps, nyears))



# Simulation

# Cell occupancy
# year 1
for (i in 1:ncells) {
  draw1 <- rmultinom(n = 1, size = 1, prob = omega[i,])
  z[i,1] <- which(draw1 == 1)
}

# years 2 and beyond
for(t in 2:nyears){
  for(i in 1:ncells){
    draw1 <- rmultinom(1,1,Phi[z[i,t-1],])
    z[i,t] <- which(draw1 == 1)
  }
}

# Local occupancy (site occupancy)
for(i in 1:ncells){
  for(j in 1:nsites){
    for(t in 1:nyears){
      draw1 <- rmultinom(1,1, theta[z[i,t],])
      u[i,j,t] <- which(draw1 == 1)
    }
  }
}

# Detection

for(i in 1:ncells){
  for(j in 1:nsites){
    for(t in 1:nyears){
      for(k in 1:nreps){
        draw1 <- rmultinom(1,1, mean.p[u[i,j,t],])
        y[i,j,k,t] <- which(draw1 == 1)
      }
    }
  }
}


##  Fit null model ------------------------------------------
# No covariates

yms <- y
str(jagsdata <- list(y = yms, ncells = dim(yms)[1], nsites = dim(yms)[2], 
                      nsurveys = dim(yms)[3],
                      nyears = dim(yms)[4]))


# Dynamic multiseason model with the parameterization for the PhiMat using growth, decline, persistence, etc

# Specify model
cat(file = "jags_txt/nulldynMSMS.txt", "
model {
  ### (1) Priors for parameters
  # State process priors
  # Priors for parameters in initial state vector (Omega)
  psi ~ dunif(0, 1)
  r ~ dunif(0, 1)
  # Priors for parameters in state transition matrix (PhiMat)
  for(s in 1:2){
    phi[s] ~ dunif(0, 1)
  }
  lgamma ~ dunif(0,1)
  G ~ dunif(0,1)
  D ~ dunif(0,1)
  
  # Priors for parameters in local occupancy (Theta)
  # theta2 ~ dunif(0, 1)              # Detection prob. when in state 2
  # for (s in 1:3) {                 # Detection prob. when in state 3
  #   beta[s] ~ dgamma(1, 1)         # Induce Dirichlet prior
  #   theta3[s] <- beta[s] / sum(beta[])
  # }
  for(s in 1:3){
    theta[s] ~ dunif(0, 1)
  }

  # Priors for parameters in observation process (varp)
  p2 ~ dunif(0, 1)                 # Detection prob. when in state 2
  for (s in 1:3) {                 # Detection prob. when in state 3
    beta[s] ~ dgamma(1, 1)         # Induce Dirichlet prior
    p3[s] <- beta[s] / sum(beta[])
  }
  ### (2) Define relationships between basic model structure and parameters
  # Define initial state vector: Year 1
  Omega[1] <- 1 - psi              # Prob. of non-occupation
  Omega[2] <- psi * (1-r)          # Prob. of occ. by a few bats
  Omega[3] <- psi * r              # Prob. of occ. by many bats
  # Define state transition probability matrix (PhiMat): years 2:nyears
  # Define probabilities of state S(t+1) given S(t)
  # For now, constant over sites and years
  # Note conditional Bernoulli parameterization of multinomial
  # Order of indices: Departing state, arrival state
  PhiMat[1,1] <- 1 - lgamma
  PhiMat[1,2] <- lgamma
  PhiMat[1,3] <- 0
  PhiMat[2,1] <- 1 - phi[1]
  PhiMat[2,2] <- phi[1] * (1 - G)
  PhiMat[2,3] <- phi[1] * G
  PhiMat[3,1] <- 1 - phi[2]
  PhiMat[3,2] <- phi[2] * D
  PhiMat[3,3] <- phi[2] * (1 - D)
  
  # Define local occupancy probability (Theta)
  Theta[1,1] <- 1
  Theta[1,2] <- 0
  Theta[1,3] <- 0
  Theta[2,1] <- 1-theta[1]
  Theta[2,2] <- theta[1]
  Theta[2,3] <- 0
  Theta[3,1] <- 1-theta[2]
  Theta[3,2] <- theta[2] * (1 - theta[3])
  Theta[3,3] <- theta[2] * theta[3]
  # Define observation probability matrix (p)
  # Order of indices: true state, observed state
  varp[1,1] <- 1
  varp[1,2] <- 0
  varp[1,3] <- 0
  varp[2,1] <- 1-p2
  varp[2,2] <- p2
  varp[2,3] <- 0
  varp[3,1] <- p3[1]
  varp[3,2] <- p3[2]
  varp[3,3] <- p3[3]
  
  ### (3) Likelihood
  # Initial state: year 1
  for (i in 1:ncells){
    z[i,1] ~ dcat(Omega[])
  }
  
  # State transitions from yearly interval
  for (i in 1:ncells){
    for (t in 2:nyears){
      z[i,t] ~ dcat(PhiMat[z[i,t-1],])
    }
  }
  
  # local occupancy
  
  for (i in 1:ncells){
    for (t in 1:nyears){
      for (j in 1:nsites){
        u[i,j,t] ~ dcat(Theta[z[i,t],])
      }
    }
  }
  
  
  # Observation equation
  for (t in 1:nyears){
    for (i in 1:ncells){
      for (j in 1:nsites){
        for (k in 1:nsurveys){
          y[i,j,k,t] ~ dcat(varp[u[i,j,t],])
        }
      }
    }
  }
  
  ### (4) Derived quantities
  # Number of sites in each state per year
  for (t in 1:nyears){
    for (i in 1:ncells){
      state1[i,t] <- equals(z[i,t], 1)   # Indicator for site in state 1
      state2[i,t] <- equals(z[i,t], 2)   # ... state 2
      state3[i,t] <- equals(z[i,t], 3)   # ... state 3
    }
    n.occ[t,1] <- sum(state1[,t])        # Number of unoccupied sites
    n.occ[t,2] <- sum(state2[,t])        # Number of sites with few bats
    n.occ[t,3] <- sum(state3[,t])        # Number of sites with many bats
    n.occ.total[t] <- n.occ[t,2] + n.occ[t, 3] # All occupied
  }
}
")

# Initial values (chosen to avoid data/model/init conflict)
zst <- array(3, dim = c(jagsdata$ncells, jagsdata$nyears) )
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("psi", "r", "phi", "lgamma", "G", "D", "theta2", "theta3", "p2", "p3", "Omega", "PhiMat",
            "Theta", "varp", "n.occ", "n.occ.total") # Could add "z"

# MCMC settings
# na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 5000 ; nc <- 3
na <- 1000 ; ni <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 3  # ~~~~ for testing, 2 mins

# Call JAGS (ART 21 min), check convergence and summarize posteriors
# odms stands for 'output dynamic multi-state'
odmsms_null <- jags(jagsdata, inits, params, "jags_txt/nulldynMSMS.txt", n.adapt = na,
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# ERROR - node inconsistent with parents
#traceplot(odms_null)
print(odmsms_null, 3)

# PhiMat
round(odmsms_null$mean$PhiMat, 2)

# Observation probabilities
round(odmsms_null$mean$Theta, 2)


# null MSMS ---------------

yms <- y[,,,1] #only one year

str(jagsdata <- list(y = yms, ncells = dim(yms)[1], nsites = dim(yms)[2], 
                      nsurveys = dim(yms)[3]))


# Dynamic multiseason model with the parameterization for the PhiMat using growth, decline, persistence, etc

# Specify model
cat(file = "jags_txt/nullMSMS.txt", "
model {
  ### (1) Priors for parameters
  # State process priors
  # Priors for parameters in initial state vector (Omega)
  psi ~ dunif(0, 1)
  r ~ dunif(0, 1)
  # # Priors for parameters in state transition matrix (PhiMat)
  # for(s in 1:2){
  #   phi[s] ~ dunif(0, 1)
  # }
  # lgamma ~ dunif(0,1)
  # G ~ dunif(0,1)
  # D ~ dunif(0,1)
  
  # Priors for parameters in local occupancy (Theta)
  for(s in 1:3){
    theta[s] ~ dunif(0, 1)
  }

  # Priors for parameters in observation process (varp)
  p2 ~ dunif(0, 1)                 # Detection prob. when in state 2
  for (s in 1:3) {                 # Detection prob. when in state 3
    beta[s] ~ dgamma(1, 1)         # Induce Dirichlet prior
    p3[s] <- beta[s] / sum(beta[])
  }
  ### (2) Define relationships between basic model structure and parameters
  # Define initial state vector: Year 1
  Omega[1] <- 1 - psi              # Prob. of non-occupation
  Omega[2] <- psi * (1-r)          # Prob. of occ. by a few bats
  Omega[3] <- psi * r              # Prob. of occ. by many bats
  # Define state transition probability matrix (PhiMat): years 2:nyears
  # Define probabilities of state S(t+1) given S(t)
  # For now, constant over sites and years
  # Note conditional Bernoulli parameterization of multinomial
  # Order of indices: Departing state, arrival state
  # PhiMat[1,1] <- 1 - lgamma
  # PhiMat[1,2] <- lgamma
  # PhiMat[1,3] <- 0
  # PhiMat[2,1] <- 1 - phi[1]
  # PhiMat[2,2] <- phi[1] * (1 - G)
  # PhiMat[2,3] <- phi[1] * G
  # PhiMat[3,1] <- 1 - phi[2]
  # PhiMat[3,2] <- phi[2] * D
  # PhiMat[3,3] <- phi[2] * (1 - D)
  
  # Define local occupancy probability (Theta)
  Theta[1,1] <- 1
  Theta[1,2] <- 0
  Theta[1,3] <- 0
  Theta[2,1] <- 1-theta[1]
  Theta[2,2] <- theta[1]
  Theta[2,3] <- 0
  Theta[3,1] <- 1-theta[2]
  Theta[3,2] <- theta[2] * (1 - theta[3])
  Theta[3,3] <- theta[2] * theta[3]
  # Define observation probability matrix (p)
  # Order of indices: true state, observed state
  varp[1,1] <- 1
  varp[1,2] <- 0
  varp[1,3] <- 0
  varp[2,1] <- 1-p2
  varp[2,2] <- p2
  varp[2,3] <- 0
  varp[3,1] <- p3[1]
  varp[3,2] <- p3[2]
  varp[3,3] <- p3[3]
  
  ### (3) Likelihood
  # Initial state: year 1
  for (i in 1:ncells){
    z[i] ~ dcat(Omega[])
  }
  
  # # State transitions from yearly interval
  # for (i in 1:ncells){
  #   for (t in 2:nyears){
  #     z[i,t] ~ dcat(PhiMat[z[i,t-1],])
  #   }
  # }
  
  # local occupancy
  
  for (i in 1:ncells){
    for (j in 1:nsites){
        u[i,j] ~ dcat(Theta[z[i],])
      }
  }
  
  
  # Observation equation
    for (i in 1:ncells){
      for (j in 1:nsites){
        for (k in 1:nsurveys){
          y[i,j,k] ~ dcat(varp[u[i,j],])
        }
      }
    }
  
  ### (4) Derived quantities
  sum.z <- sum(z[])
}
")

# Initial values (chosen to avoid data/model/init conflict)
zst <- array(3, dim = c(jagsdata$ncells) )
inits <- function(){list(z = zst)}

# Parameters monitored
params <- c("psi", "r", "theta2", "theta3", "p2", "p3", "Omega",
            "Theta", "varp", "sum.z") # Could add "z"

# MCMC settings
# na <- 1000 ; ni <- 10000 ; nt <- 5 ; nb <- 5000 ; nc <- 3
na <- 1000 ; ni <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 3  # ~~~~ for testing, 2 mins

# Call JAGS (ART 21 min), check convergence and summarize posteriors
# odms stands for 'output dynamic multi-state'
omsms_null <- jags(jagsdata, inits, params, "jags_txt/nullMSMS.txt", n.adapt = na,
                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# ERROR - node inconsistent with parents
# also the error for node inconsistent with parents










