# Analysis of Host Switching Data for Experiment 2

options(scipen = 999)

# libraries
library(tidyverse)

# import data
host_data <- read_csv("host_switch_data_averaged.csv")
cm_data_avg <- read_excel("submission1_averaged_community_matrix.xlsx")

# create data frame for total sample reads
seq_count <- cm_data_avg %>%
  filter(grepl("T23", sample)) %>%
  filter(!grepl("T23_30", sample)) %>%
  select(c(1:3, 5:6, 9, 22, 29, 32)) %>% # selecting lake trout, white sucker, and sea lamprey 
  mutate(total_reads = rowSums(across(-1)))


#------ Data Organization
# Create a detection column for each host

# add in total reads from seq_count, along with an rra column
host_data <- host_data %>%
  left_join(seq_count[, c(1, 10)], by = "sample") %>%
  mutate(
    rra_lt = round(lake_trout_reads / total_reads, 3), 
    rra_ws = round(white_sucker_reads / total_reads, 3)
  )

# Thresholds for detection:
# ----- Higher than 10 total reads
# ----- Relative read abundance higher than 1% (of total sample reads)

threshold <- 10 # 10 total read count
rra <- .01  # 1% of total reads

# add detection columns
host_data <- host_data %>%
  mutate(lt_det = ifelse(lake_trout_reads <= threshold | rra_lt < rra, 0, 1)) %>%
  mutate(ws_det = ifelse(white_sucker_reads <= threshold | rra_ws < rra, 0, 1)) %>%
  select(1:3, lt_det, ws_det, 4:9)


#------- Descriptive Stats
# can look into detections a bit more here

# Questions:
#--- How often was only the first/second host detected?
#--- How often were both hosts detected?
#--- Is there correlation between fasting length and host detection?






