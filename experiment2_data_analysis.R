# Analysis of Host Switching Data for Experiment 2

options(scipen = 999)

# libraries
library(tidyverse)
library(readxl)
library(patchwork)

# import data
host_data <- read_csv("host_switch_data_averaged.csv")
cm_data_avg <- read_excel("submission1_averaged_community_matrix.xlsx")
cm_data_reps <- read_csv("host_switch_data_replicates.csv")

# create data frame for total sample reads (average)
seq_count <- cm_data_avg %>%
  filter(grepl("T23", sample)) %>%
  filter(!grepl("T23_30", sample)) %>%
  select(c(1:3, 5:6, 9, 22, 29, 32)) %>% # selecting lake trout, white sucker, and sea lamprey 
  mutate(total_reads = rowSums(across(-1)))

# structure adjustments
host_data$fasting_days <- as.numeric(host_data$fasting_days)


#------ Data Organization (REPLICATES COMMUNITY MATRIX)

# reorder data frame
cm_data_reps <- cm_data_reps[,c(6,2:5,7:13)]
# add in total reads (both reps)
cm_data_reps <- cm_data_reps %>%
  mutate(
    total_reads_rep1 = lake_trout_reads_rep1 + white_sucker_reads_rep1 + sea_lamprey_reads_rep1,
    total_reads_rep2 = lake_trout_reads_rep2 + white_sucker_reads_rep2 + sea_lamprey_reads_rep2
  )

# add in rra counts
cm_data_reps <- cm_data_reps %>% 
  mutate(
    rra_lt_rep1 = round(lake_trout_reads_rep1 / total_reads_rep1, 3), 
    rra_ws_rep1 = round(white_sucker_reads_rep1 / total_reads_rep1, 3),
    rra_lt_rep2 = round(lake_trout_reads_rep2 / total_reads_rep2, 3), 
    rra_ws_rep2 = round(white_sucker_reads_rep2 / total_reads_rep2, 3)
)

# now need to make detection columns
# set threshold and rra
threshold <- 10 # 10 total read count
rra <- .01  # 1% of total reads

# add detection columns
host_data_reps <- cm_data_reps %>%
  mutate(lt_det_rep1 = ifelse(lake_trout_reads_rep1 <= threshold | rra_lt_rep1 < rra, 0, 1)) %>%
  mutate(ws_det_rep1 = ifelse(white_sucker_reads_rep1 <= threshold | rra_ws_rep1 < rra, 0, 1)) %>%
  mutate(lt_det_rep2 = ifelse(lake_trout_reads_rep2 <= threshold | rra_lt_rep2 < rra, 0, 1)) %>%
  mutate(ws_det_rep2 = ifelse(white_sucker_reads_rep2 <= threshold | rra_ws_rep2 < rra, 0, 1)) 




#----- Incorporate host detection data (first, second, both) into data frame

#------ REPLICATE COMMUNITY MATRIX
# set up new columns
host_data_reps <- host_data_reps %>%
  mutate(host1_det_rep1 = NA, 
         host2_det_rep1 = NA, 
         host1_det_rep2 = NA, 
         host2_det_rep2 = NA,
         both_host_det_rep1 = NA,
         both_host_det_rep2 = NA)

# loop to set detections (0 = no detection, 1 = detection)
# assign 1 or 0 to each of the host detection columns 
for (i in 1:nrow(host_data_reps)) {
  # skip the row if first or second host is NA
  if (is.na(host_data_reps$host_1[i]) & is.na(host_data_reps$host_2[i])) {
    next
  }
  
  # also skip if detection data is NA
  if (is.na(host_data_reps$lt_det_rep1[i]) & 
      is.na(host_data_reps$ws_det_rep1[i]) &
      is.na(host_data_reps$lt_det_rep2[i]) & 
      is.na(host_data_reps$ws_det_rep2[i])) {
    next
  }
  
  # add 1 if first host was detected for either replicate (and is not NA)
  if (!is.na(host_data_reps$host_1[i])) {
    # for lake trout
    if (host_data_reps$host_1[i] == "Lake Trout") {
      # check for lake trout rep 1
      if (host_data_reps$lt_det_rep1[i] == 1) {
        host_data_reps$host1_det_rep1[i] <- 1
      } else {
        host_data_reps$host1_det_rep1[i] <- 0
      }
      # check for lake trout rep 2
      if (host_data_reps$lt_det_rep2[i] == 1) {
        host_data_reps$host1_det_rep2[i] <- 1
      } else {
        host_data_reps$host1_det_rep2[i] <- 0
      }
    } 
    # and for white sucker
    if (host_data_reps$host_1[i] == "White Sucker") {
      # check for white sucker rep 1
      if (host_data_reps$ws_det_rep1[i] == 1) {
        host_data_reps$host1_det_rep1[i] <- 1
      } else {
        host_data_reps$host1_det_rep1[i] <- 0
      }
      # check for white sucker rep 2
      if (host_data_reps$ws_det_rep2[i] == 1) {
        host_data_reps$host1_det_rep2[i] <- 1
      } else {
        host_data_reps$host1_det_rep2[i] <- 0
      }
    }
  }
  
  # add 1 if second host was detected for either replicate (and is not NA)
  if (!is.na(host_data_reps$host_2[i])) {
    # for lake trout
    if (host_data_reps$host_2[i] == "Lake Trout") {
      # check for lake trout rep 1
      if (host_data_reps$lt_det_rep1[i] == 1) {
        host_data_reps$host2_det_rep1[i] <- 1
      } else {
        host_data_reps$host2_det_rep1[i] <- 0
      }
      # check for lake trout rep 2
      if (host_data_reps$lt_det_rep2[i] == 1) {
        host_data_reps$host2_det_rep2[i] <- 1
      } else {
        host_data_reps$host2_det_rep2[i] <- 0
      }
    } 
    # and for white sucker
    if (host_data_reps$host_2[i] == "White Sucker") {
      # check for white sucker rep 1
      if (host_data_reps$ws_det_rep1[i] == 1) {
        host_data_reps$host2_det_rep1[i] <- 1
      } else {
        host_data_reps$host2_det_rep1[i] <- 0
      }
      # check for white sucker rep 2
      if (host_data_reps$ws_det_rep2[i] == 1) {
        host_data_reps$host2_det_rep2[i] <- 1
      } else {
        host_data_reps$host2_det_rep2[i] <- 0
      }
    }
  }
  
  # add one if both hosts were detected (replicate 1)
  if (host_data_reps$lt_det_rep1[i] == 1 & host_data_reps$ws_det_rep1[i] == 1) {
    host_data_reps$both_host_det_rep1[i] <- 1
  } else {
    host_data_reps$both_host_det_rep1[i] <- 0
  }
  
  # add one if both hosts were detected (replicate 2)
  if (host_data_reps$lt_det_rep2[i] == 1 & host_data_reps$ws_det_rep2[i] == 1) {
    host_data_reps$both_host_det_rep2[i] <- 1
  } else {
    host_data_reps$both_host_det_rep2[i] <- 0
  }
}



#------- Descriptive Statistics

# A primary focus here is the detection of host 1 species. This indicates the
# possibility, and probabiity, that a feeding history can be detected in lamprey.
# Relating this to fasting period too (with variables such as weight gain included)
# is important too understand real-world applicability

# what were the total number and proportion of host 1 detections?
host1_rep1_det_table <- table(host_data_reps$host1_det_rep1, useNA = "no")
host1_rep2_det_table <- table(host_data_reps$host1_det_rep2, useNA = "no")

prop.table(host1_rep1_det_table)
prop.table(host1_rep2_det_table)

# looks like for both replicates, we can detect the first host about 
# 33% of the time

# comparison for the second host
host2_rep1_det_table <- table(host_data_reps$host2_det_rep1, useNA = "no")
host2_rep2_det_table <- table(host_data_reps$host2_det_rep2, useNA = "no")

prop.table(host2_rep1_det_table)
prop.table(host2_rep2_det_table)

# can detect second host about 60% of the time

# comparison for the both hosts
both_rep1_det_table <- table(host_data_reps$both_host_det_rep1, useNA = "no")
both_rep2_det_table <- table(host_data_reps$both_host_det_rep2, useNA = "no")

prop.table(both_rep1_det_table)
prop.table(both_rep2_det_table)

# for detecting both, there was less consensus, with 17% in rep 1 and 11% in rep 2



#----- Statistical Comparisons

# differences in detections between host order using logistic regression
m1_rep1 <- glm(host1_det_rep1 ~ fasting_days, data = host_data_reps, family = binomial)
m1_rep2 <- 


# Summarize the model to check the results
summary(model_1)
summary(model_2)
summary(model_both)











#-------------- OLD CODE USING AVERAGED COMMUNITY MATRIX -------------- 





#------ Data Organization (AVERAGED COMMUNITY MATRIX - OLD, USING REPLICATE DATA)
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

#------ AVERAGED COMMUNITY MATRIX (OLD - PRIMARILY USING REPLICATE DATA)
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






