# FORMATTING BAT CALL DATA AND FITTING MULTI-STATE OCCUPANCY MODELS
# Questions: f.javierarudolph@gmail.com

# Three data products were shared by Christian Stratton and the description of those is as follows:
#
# > *SiteDataClean* are data on sites including site type, clutter distance, and water distance. Also the detector type and microphone are listed. We did not get 123 data for some sites. I manually classified them based on the site description, but didn't assign a distance to water or clutter. I also included the serial numbers for the detector and microphone. Not sure if these are useful. If you are going to use them I can clean them up as there is some variability in how they were recorded between observers.
#
# > *AllCallsCapHistory* is a summary of all calls recorded by cell, site, and year and Julian day (sample night). NA values indicate that he detector did not record on that day. 0's are no calls recorded but the detector was deployed.
#
# > *WNSCallsCapHistory* is a summary of all calls that have a mean characteristic frequency at or above 34kHz.



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