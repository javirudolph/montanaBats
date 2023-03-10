# Creating this script for the skeleton of the dynamic multi-state predictions


# In this multistate model we have three different states
    # 1 - no bats
    # 2 - few bats
    # 3 - many bats

# The state of a given cell is given by the probability of it being occupied (psi)
# And if occupied, the probability of it being occupied by many bats (r)
# The vector of this probabilities is defined as Omega
# Omega is the vector of probabilities given to a multinomial distribution

# params should be given as a list and probably be the output of a different function
# So that each item in the list can be called from this function
# Also give the alternative to give no parameters and it runs a simulation

dyn_multistate_fx <- function(nyears = 5, ncells = 12, params = NULL){
  # set initial parameters
  # nyears <- 5
  # ncells <- 12
  psi <- runif(ncells, 0, 0.5)
  r <- runif(ncells, 0.3, 0.9)
  
  # State probabilities for year 1
  Omega <- matrix(c(1-psi, psi*(1-r), psi*r), ncol = 3)
  
  z <- array(NA, dim = c(ncells, nyears))
  
  # Calculate state for each cell in year 1
  for(i in 1:ncells){
    draw1 <- rmultinom(1,1,Omega[i,])
    z[i,1] <- which(draw1 == 1)
  }
  
  # Define transition probability matrix (Phi)
  # for now, just using the same matrix for all sites and all years.
  gamma1 <- runif(1)
  gamma2 <- runif(1)
  phi1 <- runif(1)
  G <- runif(1)
  phi2 <- runif(1)
  D <- runif(1)
  
  Phi <- matrix(c(1-gamma1, gamma1*(1-gamma2), gamma1*gamma2,
                  1-phi1, phi1*(1-G), phi1*G,
                  1-phi2, phi2*D, phi2*(1-D)), nrow = 3, byrow = TRUE)
  
  # fill the next years
  for(t in 2:nyears){
    for(i in 1:ncells){
      draw1 <- rmultinom(1,1,Phi[z[i,t-1],])
      z[i,t] <- which(draw1 == 1)
    }
  }
  
  return(z)
}


dyn_multistate_fx()


# This transition matrix is for site specific
gamma1 <- runif(ncells)
gamma2 <- runif(ncells)
phi1 <- runif(ncells)
G <- runif(ncells)
phi2 <- runif(ncells)
D <- runif(ncells)

# fill the next years
for(t in 2:nyears){
  for(i in 1:ncells){
    
    Phi <- matrix(c(1-gamma1[i], gamma1[i]*(1-gamma2[i]), gamma1[i]*gamma2[i],
                    1-phi1[i], phi1[i]*(1-G[i]), phi1[i]*G[i],
                    1-phi2[i], phi2[i]*D[i], phi2[i]*(1-D[i])), nrow = 3, byrow = TRUE)
    
    draw1 <- rmultinom(1,1,Phi[z[i,t-1], ])
    z[i,t] <- which(draw1 == 1)
  }
}


# What if we are changing actions yearly? 
# Then this transition matrix also changes... would have a transition matrix for each site and for each year
# Better idea... have each parameter be a matrix, where we have a parameter for each site and for each year
# Then we call the matrix for that with indices

# This transition matrix is for site and year specific
gamma1 <- matrix(runif(ncells*nyears), ncol = nyears)
gamma2 <- matrix(runif(ncells*nyears), ncol = nyears)
phi1 <- matrix(runif(ncells*nyears), ncol = nyears)
G <- matrix(runif(ncells*nyears), ncol = nyears)
phi2 <- matrix(runif(ncells*nyears), ncol = nyears)
D <- matrix(runif(ncells*nyears), ncol = nyears)


# fill the next years
for(t in 2:nyears){
  for(i in 1:ncells){
    
    Phi <- matrix(c(1-gamma1[i,t-1], gamma1[i,t-1]*(1-gamma2[i,t-1]), gamma1[i,t-1]*gamma2[i,t-1],
                    1-phi1[i,t-1], phi1[i,t-1]*(1-G[i,t-1]), phi1[i,t-1]*G[i,t-1],
                    1-phi2[i,t-1], phi2[i,t-1]*D[i,t-1], phi2[i,t-1]*(1-D[i,t-1])), nrow = 3, byrow = TRUE)
    
    draw1 <- rmultinom(1,1,Phi[z[i,t-1], ])
    z[i,t] <- which(draw1 == 1)
  }
}
