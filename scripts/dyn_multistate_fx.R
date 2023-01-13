# Creating this script for the skeleton of the dynamic multi-state predictions


# In this multistate model we have three different states
    # 1 - no bats
    # 2 - few bats
    # 3 - many bats

# The state of a given cell is given by the probability of it being occupied (psi)
# And if occupied, the probability of it being occupied by many bats (r)
# The vector of this probabilities is defined as Omega
# Omega is the vector of probabilities given to a multinomial distribution

dyn_multistate_fx <- function(nyears = 5, ncells = 12, psi = NULL, r = NULL, 
                              gamma1 = NULL, gamma2 = NULL, phi1 = NULL, G = NULL,
                              phi2 = NULL, D = NULL){
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

Phi_array <- array(NA, dim = c(3,3,ncells))
Phi_array[1,1,] <- 1-gamma1-gamma2
Phi_array[1,2,] <- gamma1
Phi_array[1,3,] <- gamma2
Phi_array[2,1,] <- 1-phi1
Phi_array[2,2,] <- phi1*(1-G)
Phi_array[2,3,] <- phi1*G
Phi_array[3,1,] <- 1-phi2
Phi_array[3,2,] <- phi2*D
Phi_array[3,3,] <- phi2*(1-D)

# fill the next years
for(t in 2:nyears){
  for(i in 1:ncells){
    draw1 <- rmultinom(1,1,Phi_array[z[i,t-1], , i])
    z[i,t] <- which(draw1 == 1)
  }
}

