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

# intercepts
alpha.lpsi <- c(0, 1.968, 1.719, 1.126, 0.150, 1.134, 0.332)
alpha.lr <- c(0, 0.115, -0.409, 0.321, 0.871, -0.465, 0.855)

# Coefficients
beta.lpsi <- c(0.64, 5.3, 1.2, 4.6, 0.15, 1.054, 0.051, -0.029)
beta.lr <- c(0.566, -0.43, 0.346)


# build the parameters psi and r based on covariates
# center and scale the covariates
center_scale <- function(x){
  center <- mean(x, na.rm = TRUE)
  x <- x - center
  scale <- sd(x, na.rm=TRUE)
  x <- x / scale
  x
}



mt_covariates %>% 
  st_drop_geometry() %>% 
  transmute(GRTS_ID,
            elev = center_scale(DEM_max),
            forest = center_scale(p_forest),
            temp = center_scale(mean_temp),
            precip = center_scale(precip),
            wetlands = center_scale(p_wetland),
            physdiv = center_scale(physio_div),
            karst = karst,
            region = as.numeric(as.factor(eco3_name))) -> scaled_covariates


scaled_covariates %>% 
  # st_drop_geometry() %>% 
  transmute(GRTS_ID,
            lpsi = alpha.lpsi[as.numeric(region)] + beta.lpsi[1] * elev + beta.lpsi[2] * elev^2 + beta.lpsi[3] * temp + beta.lpsi[4] * temp^2 + beta.lpsi[5] * physdiv + beta.lpsi[6] * precip + beta.lpsi[7] * forest + beta.lpsi[8] * wetlands,
            psi = exp(lpsi)/(1+exp(lpsi)),
            lr = alpha.lr[as.numeric(region)] + beta.lr[1] * karst + beta.lr[2] * forest + beta.lr[3] * physdiv,
            r = exp(lr)/(1+exp(lr))) -> psi_and_r

psi <- psi_and_r$psi
r <- psi_and_r$r
Omega <- matrix(c(1-psi, psi*(1-r), psi*r), ncol = 3)

ncells <- nrow(scaled_covariates)
nyears <- 5
z <- array(NA, dim = c(ncells, nyears))

# Calculate state for each cell in year 1
for(i in 1:ncells){
  draw1 <- rmultinom(1,1,Omega[i,])
  z[i,1] <- which(draw1 == 1)
}

# Viz year 1
cbind(grid_geometry, z1 = z[,1]) %>% 
  mutate(z1 = as.factor(z1)) %>% 
  ggplot() +
  geom_sf(aes(fill = z1, color = z1))




