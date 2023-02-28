#   montanaBats Project
#   Management of White Nose Syndrome


#   Code written by F. Javiera Rudolph, Ph.D.

library(tidyverse)
library(cowplot)
library(sf)



out_null_msms <- readRDS(file = "modelouts/out_null_msms.RDS")
out_M1_msms <- readRDS(file = "modelouts/out_M1_msms.RDS")
out_M2_msms <- readRDS(file = "modelouts/out_M2_msms.RDS")