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

# Need to figure out this ecoregion thing and how it is coded in the model

m1_alphas_psi <- tibble(
 region = c(),
 means = out_M1_msms$mean$alpha.lpsi,
 sds = out_M1_msms$sd$alpha.lpsi,
 model = "M1"
)


# Betas for psi
m1_betas_psi <- tibble(
  covar = c("elev", "elev_sqrd", "temp", "temp_sqrd", "physdiv", "precip", "forest", "wetlands"),
  means = out_M1_msms$mean$beta.lpsi,
  sds = out_M1_msms$sd$beta.lpsi,
  model = "M1"
  )

m2_betas_psi <- tibble(
  covar = c("elev", "elev_sqrd", "temp", "temp_sqrd", "physdiv", "precip", "forest", "wetlands"),
  means = out_M2_msms$mean$beta.lpsi,
  sds = out_M2_msms$sd$beta.lpsi,
  model = "M2"
)


betas_psi_df <- bind_rows(m1_betas_psi, m2_betas_psi) %>% 
  mutate(coeff = "psi")

betas_psi_df %>% 
  ggplot(., aes(x = means, y = covar, color = model)) +
  geom_point() +
  geom_errorbarh(aes(xmin = means - sds, xmax = means + sds)) +
  theme_bw() -> psi_fig


# Betas for r

m1_betas_r <- tibble(
  covar = c("karst", "forest", "physdiv"),
  means = out_M1_msms$mean$beta.lr,
  sds = out_M1_msms$sd$beta.lr,
  model = "M1"
)

m2_betas_r <- tibble(
  covar = c("karst", "forest", "physdiv"),
  means = out_M2_msms$mean$beta.lr,
  sds = out_M2_msms$sd$beta.lr,
  model = "M2"
)


betas_r_df <- bind_rows(m1_betas_r, m2_betas_r) %>% 
  mutate(coeff = "r")

betas_r_df %>% 
  ggplot(., aes(x = means, y = covar, color = model)) +
  geom_point() +
  geom_errorbarh(aes(xmin = means - sds, xmax = means + sds)) +
  theme_bw() -> r_fig


plot_grid(psi_fig + theme(legend.position = "none"), r_fig, rel_widths = c(0.85, 1))

model_coefficients_psi_r <- bind_rows(betas_psi_df, betas_r_df)

saveRDS(model_coefficients_psi_r, file = "modelouts/model_coefficients_psi_r.RDS")

# Make Montana map with medians for bat calls for each of the cells sampled
# compare that map to a prediction using the model for all of montana
