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

# Regions, from lines 168-170 from the estimation script:
# 1 = Northeast
# 2 = Northwest
# 3 = south


tibble(
  type = "alpha",
  covar = c("NE", "NW", "S"),
  means_m1_psi = out_M1_msms$mean$alpha.lpsi,
  sd_m1_psi = out_M1_msms$sd$alpha.lpsi,
  means_m2_psi = out_M2_msms$mean$alpha.lpsi,
  sd_m2_psi = out_M2_msms$sd$alpha.lpsi,
  means_m1_r = out_M1_msms$mean$alpha.lr,
  sds_m1_r = out_M1_msms$sd$alpha.lr,
  means_m2_r = out_M2_msms$mean$alpha.lr,
  sds_m2_r = out_M2_msms$sd$alpha.lr
) %>% 
  pivot_longer(-c(covar, type)) %>% 
  mutate(model = paste(str_match(name, "m1|m2")),
         param = paste(str_match(name, "psi|r")),
         metric = paste(str_match(name, "mean|sd"))) %>% 
  select(-name) -> alphas


tibble(
  type = "beta",
  covar = c("elev", "elev_sqrd", "temp", "temp_sqrd", "physdiv", "precip", "forest", "wetlands"),
  means_m1_psi = out_M1_msms$mean$beta.lpsi,
  sd_m1_psi = out_M1_msms$sd$beta.lpsi,
  means_m2_psi = out_M2_msms$mean$beta.lpsi,
  sd_m2_psi = out_M2_msms$sd$beta.lpsi,
) %>% 
  pivot_longer(-c(covar, type)) %>% 
  mutate(model = paste(str_match(name, "m1|m2")),
         param = paste(str_match(name, "psi|r")),
         metric = paste(str_match(name, "mean|sd"))) %>% 
  select(-name) -> betas_psi


tibble(
  type = "beta",
  covar = c("karst", "forest", "physdiv"),
  means_m1_r = out_M1_msms$mean$beta.lr,
  sds_m1_r = out_M1_msms$sd$beta.lr,
  means_m2_r = out_M2_msms$mean$beta.lr,
  sds_m2_r = out_M2_msms$sd$beta.lr
) %>% 
  pivot_longer(-c(covar, type)) %>% 
  mutate(model = paste(str_match(name, "m1|m2")),
         param = paste(str_match(name, "psi|r")),
         metric = paste(str_match(name, "mean|sd"))) %>% 
  select(-name) -> betas_r


model_coefficients_psi_r <- bind_rows(alphas, betas_psi, betas_r) %>% 
                            pivot_wider(names_from = metric, values_from = value) %>% 
                            mutate(xmin = mean - sd,
                                   xmax = mean + sd)

saveRDS(model_coefficients_psi_r, file = "modelouts/model_coefficients_psi_r.RDS")



# Visualize these coefficients to compare for the models

model_coefficients_psi_r %>% 
  filter(param == "psi") %>% 
  ggplot(., aes(x = mean, y = covar, xmin = xmin, xmax = xmax, color = model)) +
  geom_point() +
  geom_errorbarh() +
  theme_bw()


model_coefficients_psi_r %>% 
  filter(param == "r") %>% 
  ggplot(., aes(x = mean, y = covar, xmin = xmin, xmax = xmax, color = model)) +
  geom_point() +
  geom_errorbarh() +
  theme_bw()


## explore the models more:


rbind(out_null_msms$mean$n.occ,
  out_M1_msms$mean$n.occ,
  out_M2_msms$mean$n.occ) %>% 
  as_tibble() %>% 
  rename(no_bats = V1,
         few_bats = V2,
         many_bats = V3) %>% 
  rownames_to_column(var = "model") %>% 
  pivot_longer(-model, names_to = "state", values_to = "n_grids") %>% 
  mutate(prop_grids = n_grids/98) %>% 
  ggplot(aes(x = state, y = prop_grids, fill = model)) +
  geom_col(position = "dodge")




