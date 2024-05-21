# Merging Metabarcoding and Sample Data: Experiment 2 - Host-Swtiching Trials


# libraries
library(tidyverse)
library(readxl)

# import data
cm_data <- read_excel("submission1_averaged_community_matrix.xlsx")


#---------- Data Organization

# filter to only experimental 2 rows
cm_data <- cm_data %>%
  filter(grepl("T23", sample)) %>%
  filter(!grepl("T23_30", sample)) # T23_30 was from experiment 1


