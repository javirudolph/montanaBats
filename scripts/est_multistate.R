# DATA MODEL
# ESTIMATION
# SINGLE SEASON MULTISCALE OCCUPANCY


# 1. Libraries --------------------------------------------------
library(tidyverse)
library(jagsUI)
library(sf)

# 2. Load bat call data csv and modify ---------------------------

# Read the raw files - focus on WNS susceptible species only (WNSCallsCapHistory.csv)
# Data is in wide format with columns = julian date
# CSY = cell site year, where cell corresponds to the NABat Cell ID

raw_wns_calls <- read_csv("datafiles/mt_batcalls/mt_calls_csvs/WNSCallsCapHistory.csv", show_col_types = FALSE)

# What do we consider the cutoff for many bats?
manybats_cut <- 100

# Transform to long format:
raw_wns_calls %>%
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
  separate(csy, c("cell", "site", "year"), remove = FALSE) %>%
  # create binary version of detection data
  mutate(dnd = as.numeric(y > 0)) %>%
  # Create state categories, y multi state
  # Asumming a cutoff of 100 bat calls as many bats
  mutate(yms = case_when(
    y == 0 ~ 1,
    y <= manybats_cut ~ 2,
    y > manybats_cut ~ 3
  )) %>%
  # add column for number of nights of recording for each site_id
  add_count(site_id, name = "n_surveys") -> wns_calls_long

# Bundle and summarize data set




y <- wns_calls_long$dnd
sum(y)/length(y) # occupancy
str( jags.data <- list(y = y, n.cell = dim(y)[1], n.sites = dim(y)[2],
                      n.surveys = dim(y)[3], covA = covA, covB = covB, covC = covC) )

# Define model
sink("model.txt")
cat("
model {
  # Priors and model for params
  int.psi ~ dunif(0,1)         # Intercept of occupancy probability
  for(t in 1:n.sites){
    int.theta[t] ~ dunif(0,1) # Intercepts availability probability
  }
  for(t in 1:n.surveys){
    int.p[t] ~ dunif(0,1)     # Intercepts detection probability
  }
  beta.lpsi ~ dnorm(0, 0.1)    # Slopes of three covariates
  beta.ltheta ~ dnorm(0, 0.1)
  beta.lp ~ dnorm(0, 0.1)
  # 'Likelihood' (or basic model structure)
  for (i in 1:n.cell){
    # Occurrence in pond i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- logit(int.psi) + beta.lpsi * covA[i]
    for (j in 1:n.sites){
      # Occurrence in sample j
      a[i,j] ~ dbern(mu.a[i,j])
      mu.a[i,j] <- z[i] * theta[i,j]
      logit(theta[i,j]) <- logit(int.theta[j]) + beta.ltheta * covB[i,j]
      for (k in 1:n.surveys){
        # PCR detection error process in sample k
        y[i,j,k] ~ dbern(mu.y[i,j,k])
        mu.y[i,j,k] <- a[i,j] * p[i,j,k]
        logit(p[i,j,k]) <- logit(int.p[k]) + beta.lp * covC[i,j,k]
      }
    }
    tmp[i] <- step(sum(a[i,])-0.1)
  }
  # Derived quantities
  sum.z <- sum(z[])   # Total # of occupied ponds in sample
  sum.a <- sum(tmp[]) # Total # of ponds with presence in <=1 of the 5 samples
} # end model
",fill=TRUE)
sink()

# Initial values
zst <- apply(y, 1, max)        # inits for presence (z)
ast <- apply(y, c(1,2), max)   # inits for availability (a)
inits <- function() list(z = zst, a = ast, int.psi = 0.5, beta.lpsi = 0)

# Parameters monitored
params <- c("int.psi", "int.theta", "int.p", "beta.lpsi", "beta.ltheta",
            "beta.lp", "sum.z", "sum.a")

# MCMC setting
ni <- 5000   ;   nt <- 2   ;   nb <- 1000   ;   nc <- 3

# Call WinBUGS and summarize posterior
out <- jagsUI::jags(data = jags.data, inits = inits, params, "model.txt",
                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
print(out, 3)
data$sum.z








