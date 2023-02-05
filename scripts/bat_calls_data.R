#   montanaBats Project
#   Management of White Nose Syndrome


#   Code written by F. Javiera Rudolph, Ph.D.

library(tidyverse)

# 1. Data sources -----------------------------------

## 1.a. Bat call Data -------------------------------
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






