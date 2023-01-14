# This is the catch all script
# When I'm testing stuff and don't know where it goes yet

# F. Javiera Rudolph
# f.javierarudolph at gmail dot com
# or submit a Github issue

library(tidyverse)
library(sf)


mt_covariates <- readRDS("datafiles/mt_covariates.RDS")

grid_geometry <- mt_covariates %>% 
  select(GRTS_ID)

ggplot(data = grid_geometry) + geom_sf()


# what are we trying to do here?
# Create a prediction for occupancy state given a set of covariates
ncells <- nrow(mt_covariates)
nyears <- 2

psi <- 0.98
r <- 0.54

# State probabilities for year 1
Omega <- matrix(c(1-psi, psi*(1-r), psi*r), ncol = 3)

z <- array(NA, dim = c(ncells, nyears))

# Calculate state for each cell in year 1
for(i in 1:ncells){
  draw1 <- rmultinom(1,1,Omega)
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

cbind(grid_geometry, z) %>% 
  mutate(X2 = as.factor(X2)) %>% 
  ggplot() +
  geom_sf(aes(color = X2, fill = X2))


