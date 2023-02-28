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



psi_and_r <- bind_rows(c("psi", out_null_msms$mean$psi, out_null_msms$sd$psi))



c(out_null_msms$mean$psi, out_null_msms$sd$psi)

# oops, didn't track these in the models
c(out_M1_msms$mean$psi, out_M1_msms$sd$psi)
c(out_M2_msms$mean$psi, out_M2_msms$sd$psi)

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


betas_psi_df <- bind_rows(m1_betas_psi, m2_betas_psi)

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


betas_r_df <- bind_rows(m1_betas_r, m2_betas_r)

betas_r_df %>% 
  ggplot(., aes(x = means, y = covar, color = model)) +
  geom_point() +
  geom_errorbarh(aes(xmin = means - sds, xmax = means + sds)) +
  theme_bw() -> r_fig


plot_grid(psi_fig + theme(legend.position = "none"), r_fig, rel_widths = c(0.85, 1))



# Make Montana map with medians for bat calls for each of the cells sampled
# compare that map to a prediction using the model for all of montana
