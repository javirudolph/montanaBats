#   montanaBats Project
#   Management of White Nose Syndrome


#   Code written by F. Javiera Rudolph, Ph.D.


# libraries
library(tidyverse)
library(cowplot)
library(sf)


# load model outputs
out_null_msms <- readRDS(file = "modelouts/out_null_msms.RDS")
out_M1_msms <- readRDS(file = "modelouts/out_M1_msms.RDS")
out_M2_msms <- readRDS(file = "modelouts/out_M2_msms.RDS")


# Extract estimated values of psi and r for the three models and viz 

# M1 betas
m1_psi <- tibble(
  covar = c("elev", "elev_sqrd", "temp", "temp_sqrd", "physdiv", "precip", "forest", "wetlands"),
  means = out_M1_msms$mean$beta.lpsi,
  sds = out_M1_msms$sd$beta.lpsi
  )

m1_psi %>% 
  ggplot(., aes(x = means, y = covar)) +
  geom_point() +
  geom_errorbarh(aes(xmin = means - sds, xmax = means + sds))

# create figures



# Make Montana map with medians for bat calls for each of the cells sampled
# compare that map to a prediction using the model for all of montana
