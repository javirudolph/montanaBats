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
table(mt_calls$yms)/nrow(mt_calls)


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






