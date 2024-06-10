# Analysis of Host Switching Data for Experiment 2

options(scipen = 999)

# libraries
library(tidyverse)
library(readxl)
library(patchwork)
library(RColorBrewer)

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
cm_data_reps$fasting_days <- as.numeric(cm_data_reps$fasting_days)

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

# add total weight gain column
host_data_reps <- host_data_reps %>%
  mutate(total_weight_gain = weight_gain_1 + weight_gain_2)


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


# let's create a table that summarizes the totals and proportion of detections
# based on fasting days for host 1
host1_detection_table <- data.frame(
  fasting_days = unique(host_data_reps_clean$fasting_days),
  detection_total = NA,
  samples = NA,
  det_proportion = NA
)

# fill in data
for (i in 1:nrow(host1_detection_table)){
  # set up detection data frame
  rep_df <- host_data_reps_clean %>%
    select(fasting_days, host1_det_rep1, host1_det_rep2) %>%
    filter(fasting_days == host1_detection_table$fasting_days[i]) %>%
    mutate(any_host_det = ifelse(host1_det_rep1 | host1_det_rep2 == 1, 1, 0))
  # add in detection total
  host1_detection_table$detection_total[i] <- sum(rep_df$any_host_det)
  # add in sample number
  host1_detection_table$samples[i] <- nrow(rep_df)
  # add in detection proportion
  host1_detection_table$det_proportion[i] <- round(host1_detection_table$detection_total[i]/host1_detection_table$samples[i], 3)
}
# reorder
host1_detection_table <- arrange(host1_detection_table, fasting_days)
host1_detection_table

# visualize with bar chart (detection proportion)
ggplot(host1_detection_table, aes(x = as.factor(fasting_days), y = det_proportion)) +
  geom_bar(stat = "identity", fill = "#383737") +
  geom_text(aes(label = paste0("n = ", samples)), 
            vjust = -0.5, # Adjust the vertical position
            color = "black", 
            size = 3.5) + # Adjust the size of the text
  labs(x = "Fasting Days", y = "Detection Proportion", title = "Detection Proportion by Fasting Days") +
  theme_classic() 


# can also visualize total detections with bar chart
ggplot(host1_detection_table, aes(x = as.factor(fasting_days), y = detection_total)) +
  geom_bar(stat = "identity", fill = "#383737") +
  geom_text(aes(label = paste0("n = ", samples, " (", det_proportion*100, "%)")), 
            vjust = -0.5, # Adjust the vertical position
            color = "black", 
            size = 3.5) + # Adjust the size of the text
  labs(x = "Fasting Days", y = "Detection Total", title = "Host 1 Detection Total by Fasting Days") +
  theme_classic() 



#----- Statistical Comparisons

# differences in detections between host order using logistic regression
m1_rep1 <- glm(host1_det_rep1 ~ fasting_days + weight_gain_1, data = host_data_reps, family = binomial)
m1_rep2 <- glm(host1_det_rep2 ~ fasting_days + weight_gain_1, data = host_data_reps, family = binomial)

m2_rep1 <- glm(host2_det_rep1 ~ fasting_days + weight_gain_2, data = host_data_reps, family = binomial)
m2_rep2 <- glm(host2_det_rep2 ~ fasting_days + weight_gain_2, data = host_data_reps, family = binomial)

# Summarize the models to check the results
summary(m1_rep1)
summary(m1_rep2)

summary(m2_rep1) # significant relationship with fasting days
summary(m2_rep2)

# also, can remove the NA values to work with clean data
# for clean data, omit rows with NA in host1_det_rep1, host1_det_rep2, fasting_days, and rel_weight_gain
host_data_reps_clean <- host_data_reps %>%
  mutate(weight_gain_dif = weight_gain_1 - weight_gain_2,
         rel_weight_gain = weight_gain_1 / (weight_gain_1 + weight_gain_2)) %>% # adds a row to check weight difference
  filter(if_all(c(host1_det_rep1, host1_det_rep2, fasting_days, weight_gain_dif), ~ !is.na(.)))

# for host 1 clean data
m1_rep1_clean <- glm(host1_det_rep1 ~ fasting_days + weight_gain_dif, data = host_data_reps_clean, family = binomial)
m1_rep2_clean <- glm(host1_det_rep2 ~ fasting_days + weight_gain_dif, data = host_data_reps_clean, family = binomial)

summary(m1_rep1_clean)
summary(m1_rep2_clean)
# no significance




#--------- Visualizing Detections

# add column for at least one detections AND column for host species
host_data_reps_clean <- host_data_reps_clean %>%
  mutate(single_host_det = ifelse(host_data_reps_clean$host1_det_rep1 | host_data_reps_clean$host1_det_rep2 == 1, 1, 0),
         host1_species = ifelse(host_data_reps_clean$host_1 == "Lake Trout", 1, 0)) # lake trout is 1, white sucker is 0







#----- OCCUPANCY MODEL
library(unmarked)

# set up a singular, relative weight gain column to use as a covariate
# use the clean data created in previous steps

# set up detection matrix (host 1, both replicates)
detections_host1 <- data.frame(rep1 = host_data_reps_clean$host1_det_rep1,
                               rep2 = host_data_reps_clean$host1_det_rep2)


# add in site covariates
site_covs <- as.data.frame(select(host_data_reps_clean, 
                                  fasting_days, 
                                  weight_gain_dif,
                                  rel_weight_gain,
                                  host1_species)) 

# set up covariate occupancy frame object
cov_occu <- unmarkedFrameOccu(y = detections_host1, siteCovs = site_covs)
summary(cov_occu)


# set up various models to compare
occu_null <- occu(formula = ~ 1 ~1, data = cov_occu)
occu_m2 <- occu(formula = ~ weight_gain_dif + rel_weight_gain + fasting_days ~1, data = cov_occu)
occu_m3 <- occu(formula = ~ weight_gain_dif * fasting_days ~1, data = cov_occu)
occu_m4 <- occu(formula = ~ fasting_days ~1, data = cov_occu)
occu_m5 <- occu(formula = ~ weight_gain_dif ~1, data = cov_occu)
occu_m6 <- occu(formula = ~ rel_weight_gain ~1, data = cov_occu)
occu_m7 <- occu(formula = ~ weight_gain_dif + fasting_days + host1_species ~1, data = cov_occu)

# set the fit
fit <- fitList('psi(.)p(.)' = occu_null,
               'psi(.)p(weight_gain_dif + rel_weight_gain + fasting_days)' = occu_m2,
               'psi(.)p(weight_gain_dif * fasting_days)' = occu_m3,
               'psi(.)p(fasting_days)' = occu_m4,
               'psi(.)p(weight_gain_dif)' = occu_m5,
               'psi(.)p(rel_weight_gain)' = occu_m6,
               'psi(.)p(total_weight_gain + fasting_days + host1_species)' = occu_m7)
modSel(fit)
# model with fasting days and weight dif not better than null, but still close

# back transform
backTransform(occu_m6, type = "state")
backTransform(occu_m6, type = "det")


# proportional weight predictions
preds <- predict(occu_m6, type ="det", new = data.frame(rel_weight_gain = seq(-2, 2, by = 0.1)))

ggplot(data = preds, aes(x = seq(-2, 2, by = 0.1), y = Predicted)) +
  geom_smooth(stat = "smooth") +
  #geom_ribbon(data = preds, aes(ymin = preds$lower, ymax = preds$upper), alpha = 0.2)
  theme_classic()


# fasting days predictions
preds_2 <- predict(occu_m4, type ="det", new = data.frame(fasting_days = c(0:45)))

ggplot(data = preds_2, aes(x = c(0:45), y = Predicted)) +
  geom_smooth(stat = "smooth") +
  geom_ribbon(data = preds_2, aes(ymin = preds_2$lower, ymax = preds_2$upper), alpha = 0.2) +
  theme_classic()










