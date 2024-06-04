# Analysis of Host Switching Data for Experiment 2

options(scipen = 999)

# libraries
library(tidyverse)
library(readxl)

# import data
host_data <- read_csv("host_switch_data_averaged.csv")
cm_data_avg <- read_excel("submission1_averaged_community_matrix.xlsx")

# create data frame for total sample reads
seq_count <- cm_data_avg %>%
  filter(grepl("T23", sample)) %>%
  filter(!grepl("T23_30", sample)) %>%
  select(c(1:3, 5:6, 9, 22, 29, 32)) %>% # selecting lake trout, white sucker, and sea lamprey 
  mutate(total_reads = rowSums(across(-1)))

# structure adjustments
host_data$fasting_days <- as.numeric(host_data$fasting_days)

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
  select(1:3, lt_det, ws_det, 4:9) # organize order


#------- Descriptive Statistics
# can look into detections a bit more here

# Questions:
#--- How often was only the first/second host detected?
#--- How often were both hosts detected?
#--- Is there correlation between fasting length and host detection?

# How often is there a detection for first, second and both hosts?
first_host_counter <- 0
second_host_counter <- 0
both_host_counter <- 0

# loop through each row
for (i in 1:nrow(host_data)) {
  # skip the row if first or second host is NA
  if (is.na(host_data$host_1[i]) | is.na(host_data$host_2[i])) {
    next
  }
  
  # also skip if detection data is NA
  if (is.na(host_data$lt_det[i]) | is.na(host_data$ws_det[i])) {
    next
  }
  
  # check conditions for first host 
  if (host_data$host_1[i] == "Lake Trout" & host_data$lt_det[i] == 1) {
    first_host_counter <- first_host_counter + 1
  } else if (host_data$host_1[i] == "White Sucker" & host_data$ws_det[i] == 1) {
    first_host_counter <- first_host_counter + 1
  }
  
  # check conditions for second host
  if (host_data$host_2[i] == "Lake Trout" & host_data$lt_det[i] == 1) {
    second_host_counter <- second_host_counter + 1
  } else if (host_data$host_2[i] == "White Sucker" & host_data$ws_det[i] == 1) {
    second_host_counter <- second_host_counter + 1
  }
  
  # check conditions for both hosts
  if (host_data$lt_det[i] == 1 & host_data$ws_det[i] == 1) {
    both_host_counter <- both_host_counter + 1
  }
}

cat("Number of first host detections:", first_host_counter,
    "\nNumber of second host detections:", second_host_counter,
    "\nNumber of both host detections:", both_host_counter)

# comparing detections to total
total_lamprey <- host_data %>%
  filter(!is.na(lt_det) | !is.na(ws_det)) %>%
  filter(!(is.na(host_1) & is.na(host_2))) %>%
  nrow()
# 63 total lamprey samples (excluding NAs for host detections and both NAs for host data)

# looking at percentages
cat("First Host Percentage:", round((first_host_counter/total_lamprey * 100), 3),
    "\nSecond Host Percentage:", round((second_host_counter/total_lamprey * 100), 3),
    "\nBoth Host Percentage:", round((both_host_counter/total_lamprey * 100), 3))

#----- Incorporate host detection data (first, second, both) into data frame

# set up new columns
host_data <- host_data %>%
  mutate(host1_det = NA, host2_det = NA, both_host_det = NA)

# assign 1 or 0 to each of the host detection columns 
for (i in 1:nrow(host_data)) {
  # skip the row if first or second host is NA
  if (is.na(host_data$host_1[i]) & is.na(host_data$host_2[i])) {
    next
  }
  
  # also skip if detection data is NA
  if (is.na(host_data$lt_det[i]) & is.na(host_data$ws_det[i])) {
    next
  }
  
  # set detection columns to 0, if passing the two NA checks
  host_data$host1_det[i] <- 0
  host_data$host2_det[i] <- 0
  host_data$both_host_det[i] <- 0
  
  # add 1 if first host was detected (and not NA)
  if (!is.na(host_data$host_1[i]) & host_data$host_1[i] == "Lake Trout" & host_data$lt_det[i] == 1) {
    host_data$host1_det[i] <- 1
  } else if (!is.na(host_data$host_1[i]) & host_data$host_1[i] == "White Sucker" & host_data$ws_det[i] == 1) {
    host_data$host1_det[i] <- 1
  } else {
    host_data$host1_det[i] <- 0
  }
  
  # add 1 if second host was detected (and not NA)
  if (!is.na(host_data$host_2[i]) & host_data$host_2[i] == "Lake Trout" & host_data$lt_det[i] == 1) {
    host_data$host2_det[i] <- 1
  } else if (!is.na(host_data$host_2[i]) & host_data$host_2[i] == "White Sucker" & host_data$ws_det[i] == 1) {
    host_data$host2_det[i] <- 1
  } else {
    host_data$host2_det[i] <- 0
  }
  
  # add one of both hosts were detected
  if (host_data$lt_det[i] == 1 & host_data$ws_det[i] == 1) {
    host_data$both_host_det[i] <- 1
  } else {
    host_data$both_host_det[i] <- 0
  }
}

#----- Looking at correlations between detections and fasting period

# remove rows with NA values in relevant columns
clean_data <- host_data %>%
  filter(!is.na(fasting_days) & !is.na(host1_det) & !is.na(host2_det) & !is.na(both_host_det))

# Calculate point biserial correlation for each binary column
cor_host1_det <- cor.test(clean_data$fasting_days, clean_data$host1_det, method = "pearson")
cor_host2_det <- cor.test(clean_data$fasting_days, clean_data$host2_det, method = "pearson")
cor_both_host_det <- cor.test(clean_data$fasting_days, clean_data$both_host_det, method = "pearson")

# Create a boxplot for host1_det
ggplot(clean_data, aes(x = factor(host1_det), y = fasting_days)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, color = "blue") +
  labs(title = "Fasting Period vs Host1 Detection", x = "Host1 Detected (0 = No, 1 = Yes)", y = "Fasting Period")

# Create a boxplot for host2_det
ggplot(clean_data, aes(x = factor(host2_det), y = fasting_days)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, color = "blue") +
  labs(title = "Fasting Period vs Host2 Detection", x = "Host2 Detected (0 = No, 1 = Yes)", y = "Fasting Period")

# Create a boxplot for both_host_det
ggplot(clean_data, aes(x = factor(both_host_det), y = fasting_days)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, color = "blue") +
  labs(title = "Fasting Period vs Both Hosts Detection", x = "Both Hosts Detected (0 = No, 1 = Yes)", y = "Fasting Period")













