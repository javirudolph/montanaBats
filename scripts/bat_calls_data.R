#   montanaBats Project
#   Management of White Nose Syndrome


#   Code written by F. Javiera Rudolph, Ph.D.

library(tidyverse)
library(cowplot)

# 1. Data sources -----------------------------------

## 1a. Bat call Data -------------------------------
# Data here is not specific to any bat species, but it is a subset of calls
# the frequency used to subset corresponds to bat species known to be 
# susceptible to white nose syndrome

raw_mtcalls <- read_csv("datafiles/mt_batcalls/mt_calls_csvs/WNSCallsCapHistory.csv")

raw_mtcalls %>% 
  # create variable for julian date, formerly in the columns
  pivot_longer(., 4:last_col(), names_to = "jdate", values_to = "y") %>%
  # keep only needed variables
  select(SiteID, CSY, jdate, y) %>%
  # julian date should be a number
  mutate(jdate = as.double(jdate)) %>%
  # remove NAs - means recorders wasn't on
  drop_na(y) %>%
  # clean the names
  janitor::clean_names() %>%
  # separate CSY into three
  separate(csy, c("cell", "site", "year")) %>%
  # create binary version of detection data
  mutate(dnd = as.numeric(y > 0)) -> mt_calls

# Summary of the bat call data

median_calls <- median(mt_calls$y)
ggplot(data = mt_calls, aes(x = y)) + 
  geom_histogram() + 
  geom_vline(xintercept = median_calls, color = "red") -> call_dist_plot

mt_calls %>% 
  mutate(yms = case_when(
    y == 0 ~ 1,
    y <= median_calls ~ 2, 
    y > median_calls ~ 3
  )) -> mt_calls

ggplot(data = mt_calls, aes(x = yms, fill = year)) +
  geom_bar() -> state_dist_plot

plot_grid(call_dist_plot, state_dist_plot)

# naive occupancy probability not accounting for imperfect detection or multilevel occupancy
# nested sites

sum(mt_calls$dnd)/nrow(mt_calls)
table(mt_calls$yms)
table(mt_calls$yms)/nrow(mt_calls)


# Number of cells, sites per cell, surveys per site:
cell_list <- sort(unique(mt_calls$cell))
ncells <- length(cell_list)

# distribution of number of surveys per site
mt_calls %>% 
  count(cell, site) %>% 
  ggplot(., aes(x = n)) +
  geom_histogram()

# How many sites per cell per year
mt_calls %>% 
  select(site_id, cell, year) %>% 
  distinct() %>% 
  count(cell, year) %>% 
  count(n)

nsites <- mt_calls %>% select(site_id, cell) %>% distinct() %>% count(cell) %>% pull(n) %>% max()

# how many surveys/nights per site:
mt_calls %>% 
  count(site_id, cell) %>% 
  pull(n) %>% 
  max() -> max_surveys

# How many cells have info for both years?
mt_calls %>% 
  select(site_id, cell, year) %>% 
  distinct() %>% 
  select(cell, year) %>% 
  distinct() %>% 
  count(cell) %>% 
  filter(n ==2) %>% 
  pull(cell) -> cell_list_2years

length(cell_list_2years)


## 1b. Site covariates -----------------

raw_sitecovs <- read_csv("datafiles/mt_batcalls/mt_calls_csvs/SiteDataCleanFinal.csv") %>% 
                    mutate(Site_Type = str_remove(Site_Type, " "))

# interested in site type as a covariate for local availability
unique(raw_sitecovs$Site_Type)

# What falls under the "other" category
raw_sitecovs %>% 
  filter(Site_Type == "Other") %>% 
  View()

raw_sitecovs %>% 
  select(SiteID, Site_Type) %>% 
  group_by(Site_Type) %>% 
  tally()

sitecov


# Models ----------------
## Null ----------------

# Multi-state multi-scale occupancy model
# no covariates

# Cell level occupancy
# z_i ~ multinom(Omega)
# Omega = [1-psi, psi(1-r), psi*r]

# Site level availability
# u_{i,j} ~ multinom(Theta)
# Theta = [1, 0, 0
#          1-theta2, theta2, 0
#          1-theta32-theta33, theta32, theta33]

# Observation process
# y_{i,j,k} ~ multinom(p)
# p = [1,0,0
#      1-p2, p2, 0,
#      1-p32-p33, p32, p33]

# Create new site identifier that pools data from both years
# and then creates an identifier for cell_site
mt_calls %>% 
  group_split(cell) %>% 
  purrr::map_df(~.x %>% group_by(site_id) %>% mutate(id = cur_group_id())) %>% 
  mutate(cell_site = paste(cell, id, sep = "_")) -> null_msms_dat


# bundle data
yms <- array(NA, dim = c(ncells, nsites, max_surveys),
             dimnames = list(cell_list, 1:nsites, 1:max_surveys))

for(i in 1:ncells){
  for(j in 1:nsites){
    sel_cell_site <- paste(cell_list[i], j, sep = "_")
    tmp <- null_msms_dat[null_msms_dat$cell_site == sel_cell_site,]
    nr <- nrow(tmp)
    if(nr > 0){
      yms[i,j, 1:nr] <- tmp$yms
    }
  }
}

# There is variation in the number of nights that microphones are recording for each cell or site
nsurveys <- array(NA, dim = c(ncells, nsites))
for(i in 1:ncells){
  for(j in 1:nsites){
    tmp <- which(!is.na(yms[i,j,]))
    if(length(tmp) > 0){
      nsurveys[i,j] <- max(tmp)
    }
  }
}

str(batdata <- list(y = yms, ncells = dim(yms)[1], nsites = dim(yms)[2], nsurveys = nsurveys))

# Specify the model

cat(file = "jags_txt/null_msms.txt", "
model {
  # Priors
  psi ~ dunif(0, 1)
  r ~ dunif(0, 1)
  theta2 ~ dunif(0,1)
  p2 ~ dunif(0, 1)
  
  # Multinomial logit link for availability model for state 3 (= many bats)
  ltheta32 ~ dnorm(0, 0.001)
  ltheta33 ~ dnorm(0, 0.001)
  theta32 <- exp(ltheta32) / (1 + exp(ltheta32) + exp(ltheta33))
  theta33 <- exp(ltheta33) / (1 + exp(ltheta32) + exp(ltheta33))
  theta31 <- 1-theta32-theta33  
  
  # Multinomial logit link for observation model for state 3
  lp32 ~ dnorm(0, 0.001)
  lp33 ~ dnorm(0, 0.001)
  p32 <- exp(lp32) / (1 + exp(lp32) + exp(lp33))
  p33 <- exp(lp33) / (1 + exp(lp32) + exp(lp33))
  p31 <- 1-p32-p33                     
  
  # Define initial state vector (Omega)
  Omega[1] <- 1 - psi                  # Prob. of no bats
  Omega[2] <- psi * (1-r)              # Prob. of occupancy (w/ few bats)
  Omega[3] <- psi * r                  # Prob. of occupancy (with many bats)
  
  # Define availability matrix (Theta)
  # Order of indices: true state, observed state
  Theta[1,1] <- 1
  Theta[1,2] <- 0
  Theta[1,3] <- 0
  Theta[2,1] <- 1-theta2
  Theta[2,2] <- theta2
  Theta[2,3] <- 0
  Theta[3,1] <- theta31                    
  Theta[3,2] <- theta32
  Theta[3,3] <- theta33
  
  # Define observation matrix (pDet)
  # Order of indices: true state, observed state
  pDet[1,1] <- 1
  pDet[1,2] <- 0
  pDet[1,3] <- 0
  pDet[2,1] <- 1-p2
  pDet[2,2] <- p2
  pDet[2,3] <- 0
  pDet[3,1] <- p31                    # = 1-p32-p33 as per prior section
  pDet[3,2] <- p32
  pDet[3,3] <- p33
  
  # State-space likelihood
  # State equation: model of true states (z)
  for (i in 1:ncells){
    z[i] ~ dcat(Omega[])
  }
  
  # Availability equation
  for (i in 1:ncells){
    for (j in 1:nsites){
      a[i,j] ~ dcat(Theta[z[i],])
    }
  }
  
  # Observation equation
  for (i in 1:ncells){
    for (j in 1:nsites){
      for (k in 1:nsurveys[i,j]){
        y[i,j,k] ~ dcat(pDet[a[i,j],])
      }
    }
  }
  
  # Derived quantities
  for (i in 1:ncells){
    occ1[i] <- equals(z[i], 1)
    occ2[i] <- equals(z[i], 2)
    occ3[i] <- equals(z[i], 3)
  }
  n.occ[1] <- sum(occ1[]) # Grid cells in state 1
  n.occ[2] <- sum(occ2[]) # Grid cells in state 2
  n.occ[3] <- sum(occ3[]) # Grid cells in state 3
}
")

# Initial values
zst <- rep(3, nrow(batdata$y)) # Initialize at highest possible state
inits <- function(){list(z = zst)}

# Parameters monitored (could add "z")
params <- c("psi", "r", "p2", "p31", "p32", "p33", "Omega", "Theta", "pDet", "n.occ")

# MCMC settings
na <- 1000 ; ni <- 2000 ; nt <- 2 ; nb <- 1000 ; nc <- 3

# Call JAGS, check convergence and summarize posteriors
library(jagsUI)
out_null_msms <- jags(batdata, inits, params, "jags_txt/null_msms.txt", n.adapt = na,
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
# par(mfrow = c(3,3))  # ~~~ no longer needed
traceplot(out1)
print(out1, 3)





