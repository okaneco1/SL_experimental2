# Analysis of Host Switching Data for Experiment 2

options(scipen = 999)

# libraries
library(tidyverse)
library(readxl)
library(patchwork)
library(RColorBrewer)
library(gridExtra)

# import data
host_data <- read_csv("host_switch_data_averaged.csv")
cm_data_avg <- read_excel("submission1_averaged_community_matrix.xlsx")
cm_data_reps <- read_csv("host_switch_data_replicates.csv")

# create data frame for total sample reads (average)
seq_count <- cm_data_avg %>%
  filter(grepl("T23", sample)) %>%
  # specific samples to remove
  filter(!grepl("T23_30", sample) & !grepl("T23_18", sample)) %>%
  select(c(1:3, 5:6, 9, 22, 29, 32)) %>% # selecting lake trout, white sucker, and sea lamprey 
  mutate(total_reads = rowSums(across(-1)))

# structure adjustments
host_data$fasting_days <- as.numeric(host_data$fasting_days)
cm_data_reps$fasting_days <- as.numeric(cm_data_reps$fasting_days)

#------ Data Organization (REPLICATES COMMUNITY MATRIX)
# filter out samples
cm_data_reps <- cm_data_reps %>%
  filter(!grepl("T23_30", sample) & !grepl("T23_18", sample))

# reorder data frame
cm_data_reps <- cm_data_reps[,c(8,2:7,9:15)]
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

# if weight loss, register as 0 (no weight gain)
host_data_reps$weight_gain_1[host_data_reps$weight_gain_1 < 0] <- 0
host_data_reps$weight_gain_2[host_data_reps$weight_gain_2 < 0] <- 0

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
         both_host_det_rep2 = NA,
         both_hosts = NA)

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
  filter(if_all(c(host1_det_rep1, host1_det_rep2, fasting_days, weight_gain_dif, days_attached_1, days_attached_2), ~ !is.na(.)))

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
         host1_species = ifelse(host_data_reps_clean$host_1 == "Lake Trout", 1, 0)) %>% # lake trout is 1, white sucker is 0
  arrange(sample)







#--------------- OCCUPANCY MODEL
library(unmarked)

# some visualizations first
long_det_data <- host_data_reps_clean %>%
  pivot_longer(cols = c(host1_det_rep1, host1_det_rep2),
               names_to = "replicate",
               values_to = "host1_detection")

# relative weight gain
ggplot(long_det_data, aes(x = rel_weight_gain, y = host1_detection, color = replicate)) +
  geom_jitter(width = 0.01, height = 0.01) +
  labs(x = "Relative Weight Gain", y = "Host 1 Detection", 
       title = "Host 1 Detections per Replicate vs Relative Weight Gain",
       color = "Replicate") +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("#87CEEB", "#228B22"),
                     labels = c("Replicate 1", "Replicate 2")) +
  geom_hline(yintercept = 0, color = "#5c5e5c", linewidth = 0.5) +   
  geom_vline(xintercept = 0, color = "#5c5e5c", linewidth = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.position = "inside",
        legend.title = element_text(face = "bold"),
        legend.position.inside = c(0.85, 0.7),
        legend.background = element_rect(fill = alpha('white', 0.5), color = NA))

# fasting days
ggplot(long_det_data, aes(x = fasting_days, y = host1_detection, color = replicate)) +
  geom_jitter(width = 0.5, height = 0.02) +
  labs(x = "Fasting Days", y = "Host 1 Detection", 
       title = "Host 1 Detections per Replicate vs Fasting Days",
       color = "Replicate") +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("#87CEEB", "#228B22"),
                     labels = c("Replicate 1", "Replicate 2")) +
  geom_hline(yintercept = 0, color = "#5c5e5c", linewidth = 0.5) +   
  geom_vline(xintercept = 0, color = "#5c5e5c", linewidth = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.position = "inside",
        legend.title = element_text(face = "bold"),
        legend.position.inside = c(0.85, 0.7),
        legend.background = element_rect(fill = alpha('white', 0.5), color = NA))



# set up a singular, relative weight gain column to use as a covariate
# use the clean data created in previous steps

# set up detection matrix (host 1, both replicates)
detections_host1 <- data.frame(rep1 = host_data_reps_clean$host1_det_rep1,
                               rep2 = host_data_reps_clean$host1_det_rep2)


# add in site covariates
site_covs <- as.data.frame(select(host_data_reps_clean, 
                                  fasting_days,
                                  rel_weight_gain,
                                  host1_species,
                                  weight_gain_2,
                                  days_attached_1,
                                  days_attached_2)) 

# set up covariate occupancy frame object
cov_occu <- unmarkedFrameOccu(y = detections_host1, siteCovs = site_covs)
summary(cov_occu)


# set up various models to compare
occu_null <- occu(formula = ~ 1 ~1, data = cov_occu)
occu_m2 <- occu(formula = ~ rel_weight_gain + fasting_days + host1_species ~1, data = cov_occu)
occu_m3 <- occu(formula = ~ rel_weight_gain * fasting_days ~1, data = cov_occu)
occu_m4 <- occu(formula = ~ fasting_days ~1, data = cov_occu)
occu_m5 <- occu(formula = ~ rel_weight_gain ~1, data = cov_occu)
occu_m6 <- occu(formula = ~ host1_species ~1, data = cov_occu)
occu_m7 <- occu(formula = ~ rel_weight_gain + fasting_days ~1, data = cov_occu)
occu_m8 <- occu(formula = ~ weight_gain_2 ~1, data = cov_occu)
occu_m9 <- occu(formula = ~ weight_gain_2 + rel_weight_gain ~1, data = cov_occu)
occu_m10 <- occu(formula = ~ weight_gain_2 + rel_weight_gain + fasting_days ~1, data = cov_occu)
occu_m11 <- occu(formula = ~ days_attached_1 ~1, data = cov_occu)
occu_m12 <- occu(formula = ~ days_attached_2 ~1, data = cov_occu)
occu_m13 <- occu(formula = ~ days_attached_2 + rel_weight_gain ~1, data = cov_occu)

# set the fit
fit <- fitList('psi(.)p(.)' = occu_null,
               'psi(.)p(rel_weight_gain + fasting_days + host1_species)' = occu_m2,
               'psi(.)p(rel_weight_gain * fasting_days)' = occu_m3,
               'psi(.)p(fasting_days)' = occu_m4,
               'psi(.)p(rel_weight_gain)' = occu_m5,
               'psi(.)p(host1_species)' = occu_m6,
               'psi(.)p(rel_weight_gain + fasting_days)' = occu_m7,
               'psi(.)p(weight_gain_2)' = occu_m8,
               'psi(.)p(weight_gain_2 + rel_weight_gain)' = occu_m9,
               'psi(.)p(weight_gain_2 + rel_weight_gain + fasting_days)' = occu_m10,
               'psi(.)p(days_attached_1)' = occu_m11,
               'psi(.)p(days_attached_2)' = occu_m12,
               'psi(.)p(days_attached_2 + rel_weight_gain)' = occu_m13)
modSel(fit)
# model with fasting days and weight dif not better than null, but still close

# back transform
#backTransform(occu_m5, type = "state")
backTransform(occu_null, type = "det") # null shows detection 

# setting up as a linear combination
det_est <- coef(occu_m7, type = "det")
# Create a linear combination for the detection estimates
# Assuming your detection formula has an intercept and two covariates (rel_weight_gain and fasting_days)
lc_det <- linearComb(occu_m7, coefficients = c(1, mean(cov_occu@siteCovs$rel_weight_gain), mean(cov_occu@siteCovs$fasting_days)), type = "det")
# Back-transform the linear combination
back_transform_det <- backTransform(lc_det)

# can also use an anti-logit function
antilogit <- function(x) { exp(x) / (1 + exp(x) ) }
antilogit(coef(occu_m5)["p(rel_weight_gain)"])


# PREDICTIONS

# proportional weight predictions
preds <- predict(occu_m5, type ="det", new = data.frame(rel_weight_gain = seq(0, 1, by = 0.1)))

ggplot(data = preds, aes(x = seq(0, 1, by = 0.1), y = Predicted)) +
  geom_smooth(stat = "smooth") +
  geom_ribbon(data = preds, aes(ymin = lower, ymax = upper), alpha = 0.2)
  theme_classic()


# fasting days predictions (with model that incorporates)
preds_2 <- predict(occu_m7, type ="det", new = data.frame(fasting_days = c(0:45),rel_weight_gain = seq(0, 1, length = 46)))

ggplot(data = preds_2, aes(x = c(0:45), y = Predicted)) +
  geom_smooth(stat = "smooth") +
  geom_ribbon(data = preds_2, aes(ymin = preds_2$lower, ymax = preds_2$upper), alpha = 0.2) +
  theme_classic()








# ------- LINEAR MODELS COMPARING READ COUNT

# adding a host 1 read count column
host_data_reps_clean$host1_read_count_rep1 <- NA
host_data_reps_clean$host1_read_count_rep2 <- NA

for (i in 1:nrow(host_data_reps_clean)) {
  if (host_data_reps_clean$host_1[i] == "Lake Trout") {
    host_data_reps_clean$host1_read_count_rep1[i] <- host_data_reps_clean$lake_trout_reads_rep1[i]
    host_data_reps_clean$host1_read_count_rep2[i] <- host_data_reps_clean$lake_trout_reads_rep2[i]
  }
  else if (host_data_reps_clean$host_1[i] == "White Sucker") {
    host_data_reps_clean$host1_read_count_rep1[i] <- host_data_reps_clean$white_sucker_reads_rep1[i]
    host_data_reps_clean$host1_read_count_rep2[i] <- host_data_reps_clean$white_sucker_reads_rep2[i]
  }
}

# viewing data
hist(host_data_reps_clean$host1_read_count_rep1)
hist(host_data_reps_clean$host1_read_count_rep2)
# reshape for plot
long_data <- host_data_reps_clean %>%
  pivot_longer(cols = c(host1_read_count_rep1, host1_read_count_rep2),
               names_to = "replicate",
               values_to = "read_count")

# plot each variable against read count (for each rep)
ggplot(long_data, aes(x = rel_weight_gain, y = read_count, color = replicate)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal()
ggplot(long_data, aes(x = fasting_days, y = read_count, color = replicate)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  geom_jitter() +
  theme_minimal()
ggplot(long_data, aes(x = host1_species, y = read_count, color = replicate)) +
  geom_point() +
  geom_smooth(method = "loess") +
  theme_minimal()


# setting up linear models
log_lm1 <- lm(log1p(host1_read_count_rep1) ~ rel_weight_gain + host1_species + fasting_days, data = host_data_reps_clean)


# checking plot
par(mfrow = c(2, 2))
plot(log_lm1)

par(mfrow = c(1, 1))

summary(log_lm1)






# ------- misc._data 
long_det_host2_data <- host_data_reps_clean %>%
  pivot_longer(cols = c(host2_det_rep1, host2_det_rep2),
               names_to = "replicate",
               values_to = "host2_detection")
# fasting days
host2_fast <- ggplot(long_det_host2_data, aes(x = fasting_days, y = host2_detection, color = replicate)) +
  geom_point() +
  geom_jitter(width = 0.5, height = 0.02) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal()
# relative weight gain (host 1)
ggplot(long_det_host2_data, aes(x = rel_weight_gain, y = host2_detection, color = replicate)) +
  geom_point() +
  geom_jitter(width = 0.01, height = 0.02) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal()

# plot both fasting graphs
grid.arrange(host1_fast, host2_fast, ncol=2)

#check
which(host_data_reps_clean$host1_det_rep1+host_data_reps_clean$host1_det_rep2 > 0 & host_data_reps_clean$host2_det_rep1+host_data_reps_clean$host2_det_rep2 > 0)
view(host_data_reps_clean[which(host_data_reps_clean$host1_det_rep1+host_data_reps_clean$host1_det_rep2 > 0 & host_data_reps_clean$host2_det_rep1+host_data_reps_clean$host2_det_rep2 > 0), ])






# ------------------------------------------------------------------
# looking at removing if hosts did NOT feed on a second host

# additionally, using pcount() to set occupancy 1 (essentially does not
# include an occupancy calculation in the model)
# ------------------------------------------------------------------

# this is because if there is no second host feeding, detection of 
# an feeding history is nullified, as only a single host was fed on

# start by simplifying the data (only lamprey with two hosts)
two_hosts_data <- host_data_reps_clean %>%
  filter(weight_gain_2 > 0)

# some visualizations first (as before)
long_two_hosts_data <- two_hosts_data %>%
  pivot_longer(cols = c(host1_det_rep1, host1_det_rep2),
               names_to = "replicate",
               values_to = "host1_detection")

# relative weight gain (host 1)
host1_two_hosts_weight_plot <- ggplot(long_two_hosts_data, aes(x = rel_weight_gain, y = host1_detection, color = replicate)) +
  geom_jitter(width = 0.01, height = 0.01) +
  labs(x = "Relative Weight Gain", y = "Host 1 Detection", 
       title = "Host 1 Detections per Replicate vs Relative Weight Gain",
       color = "Replicate") +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("#87CEEB", "#228B22"),
                     labels = c("Replicate 1", "Replicate 2")) +
  geom_hline(yintercept = 0, color = "#5c5e5c", linewidth = 0.5) +   
  geom_vline(xintercept = 0, color = "#5c5e5c", linewidth = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.position = "inside",
        legend.title = element_text(face = "bold"),
        legend.position.inside = c(0.85, 0.7),
        legend.background = element_rect(fill = alpha('white', 0.5), color = NA))


# optional to remove LOESS warnings, adds jitter to fasting days so they are
# treated more as a continuous variable

#long_two_hosts_data$fasting_days <- jitter(long_two_hosts_data$fasting_days, factor=1)

# fasting days
host1_two_hosts_fasting_plot <- ggplot(long_two_hosts_data, aes(x = fasting_days, y = host1_detection, color = replicate)) +
  geom_jitter(width = 0.5, height = 0.02) +
  labs(x = "Fasting Days", y = "Host 1 Detection", 
       title = "Host 1 Detections per Replicate",
       color = "Replicate") +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("#87CEEB", "#228B22"),
                     labels = c("Replicate 1", "Replicate 2")) +
  geom_hline(yintercept = 0, color = "#5c5e5c", linewidth = 0.5) +   
  geom_vline(xintercept = 0, color = "#5c5e5c", linewidth = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.position = "inside",
        legend.title = element_text(face = "bold"),
        legend.position.inside = c(0.2, 0.7),
        legend.background = element_rect(fill = alpha('white', 0.5), color = NA))

# NOTE: warnings do not seem to make a difference compared to non-jittered graph


# set up detection matrix (host 1, both replicates, both host data)
detections_host1_two_hosts <- data.frame(rep1 = two_hosts_data$host1_det_rep1,
                                         rep2 = two_hosts_data$host1_det_rep2)

# add in site covariates (simplified from previous)
site_covs2 <- as.data.frame(select(two_hosts_data, 
                                  fasting_days,
                                  rel_weight_gain)) 

# set up covariate occupancy frame object
cov_occu2 <- unmarkedFramePCount(y = detections_host1_two_hosts, siteCovs = site_covs2)
summary(cov_occu2)


# set up various models to compare (just with relative weight gain and fasting period)
pcount_null_2 <- pcount(~ 1 ~ rel_weight_gain + fasting_days, data = cov_occu2, K = 3)
pcount_m1_2 <- pcount(~ rel_weight_gain ~ rel_weight_gain + fasting_days, data = cov_occu2, K = 3)
pcount_m2_2 <- pcount(~ rel_weight_gain + fasting_days ~ rel_weight_gain + fasting_days, data = cov_occu2, K = 3)

# set the fit
fit2 <- fitList('p(.)' = pcount_null_2,
                'p(rel_weight_gain)' = pcount_m1_2,
                'p(rel_weight_gain + fasting_days)' = pcount_m2_2)
modSel(fit2)
# model with fasting days and weight dif not better than null, but still close


# PREDICTIONS
# proportional weight predictions
preds <- predict(pcount_m2_2, type ="det", newdata = data.frame(rel_weight_gain = seq(0, 1, length = 46),
                                                                fasting_days = seq(0, 45, by =1)))

# plot predictions
ggplot(data = preds, aes(x = seq(0, 1, length = 46), y = Predicted)) +
  geom_smooth(stat = "smooth") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  theme_classic()


# -----
# contour plot of occupancy model
library(geomtextpath)

# generate grid points
ngrid <- 45   # Number of grid points
xvals <- seq(from = min(site_covs$fasting_days), to = max(site_covs$fasting_days), length.out = ngrid + 1)
yvals <- seq(from = min(site_covs$rel_weight_gain), to = max(site_covs$rel_weight_gain), length.out = ngrid + 1)

# generate matrix of occupancy probabilities
x.elev <- rep(xvals, length(xvals))
y.length <- rep(yvals, each = length(xvals))

# prediction function
zvals <- matrix(predict(pcount_m2_2, type="det", new=data.frame(fasting_days=x.elev, rel_weight_gain=y.length))$Predict, nrow=ngrid+1)

# make data suitable for ggplot
contour_data <- expand.grid(fasting_days = xvals, rel_weight_gain = yvals)
contour_data$occupancy_prob <- as.vector(zvals)

# contour plot
ggplot(contour_data, aes(x = fasting_days, y = rel_weight_gain, z = occupancy_prob)) +
  stat_contour_filled(breaks = seq(0, 1, by = 0.1)) +
  geom_textcontour(linecolour = "#242424", textcolour = "#242424") + 
  labs(x = "Fasting Days", y = "Relative Weight Gain (Host 1)", 
       title = "Host 1 Predicted Detection Probabilities") +  
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12))





#---------------------------
# visualizing host 1 and host 2 together for two hosts data

long_two_hosts_data_2 <- host_data_reps_clean %>%
  filter(weight_gain_2 > 0) %>%
  pivot_longer(cols = c(host2_det_rep1, host2_det_rep2),
               names_to = "replicate",
               values_to = "host2_detection")

# relative weight gain (host 2)
host2_two_hosts_weight_plot <- ggplot(long_two_hosts_data_2, aes(x = rel_weight_gain, y = host2_detection, color = replicate)) +
  geom_jitter(width = 0.01, height = 0.01) +
  labs(x = "Relative Host 1 Weight Gain", y = "Host 2 Detection", 
       title = "Host 2 Detections per Replicate vs Relative Host 1 Weight Gain",
       color = "Replicate") +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("#87CEEB", "#228B22"),
                     labels = c("Replicate 1", "Replicate 2")) +
  geom_hline(yintercept = 0, color = "#5c5e5c", linewidth = 0.5) +   
  geom_vline(xintercept = 0, color = "#5c5e5c", linewidth = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.position = "inside",
        legend.title = element_text(face = "bold"),
        legend.position.inside = c(0.85, 0.7),
        legend.background = element_rect(fill = alpha('white', 0.5), color = NA))

# fasting days (host 2)
host2_two_hosts_fasting_plot <- ggplot(long_two_hosts_data_2, aes(x = fasting_days, y = host2_detection, color = replicate)) +
  geom_jitter(width = 0.5, height = 0.02) +
  labs(x = "Fasting Days", y = "Host 2 Detection", 
       title = "Host 2 Detections per Replicate",
       color = "Replicate") +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("#87CEEB", "#228B22"),
                     labels = c("Replicate 1", "Replicate 2")) +
  geom_hline(yintercept = 0, color = "#5c5e5c", linewidth = 0.5) +   
  geom_vline(xintercept = 0, color = "#5c5e5c", linewidth = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.position = "inside",
        legend.title = element_text(face = "bold"),
        legend.position.inside = c(0.2, 0.35),
        legend.background = element_rect(fill = alpha('white', 0.5), color = NA))


# plot both fasting graphs
grid.arrange(host1_two_hosts_fasting_plot, host2_two_hosts_fasting_plot, ncol=2)

# can create additional column for both detections (in any replicate)
# finally, add one for both host detections of either replicate
two_hosts_data$both_hosts <- 0

# count both host detections
for (i in 1:nrow(two_hosts_data)) {
  if (two_hosts_data$host1_det_rep1[i] == 1 | two_hosts_data$host1_det_rep2[i] == 1) {
    if (two_hosts_data$host2_det_rep1[i] == 1 | two_hosts_data$host2_det_rep2[i] == 1) {
      two_hosts_data$both_hosts[i] <- 1
    }
  }
}

# fasting plot for both hosts detections
both_hosts_fasting_plot <- ggplot(host_data_reps_clean, aes(x = fasting_days, y = both_hosts)) +
  geom_jitter(width = 0.5, height = 0.02, color = "#003582") +
  labs(x = "Fasting Days", y = "Both Hosts Detection", 
       title = "Both Host Detections in Either Replicate") +
  geom_smooth(method = "loess", se = FALSE, color = "#003582") +
  geom_hline(yintercept = 0, color = "#5c5e5c", linewidth = 0.5) +   
  geom_vline(xintercept = 0, color = "#5c5e5c", linewidth = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12))

# plot all together
grid.arrange(host1_two_hosts_fasting_plot, 
             host2_two_hosts_fasting_plot,
             both_hosts_fasting_plot, ncol=2)



#----------------------------------------------------------------------------
# looking at occupancy models and contour plots for both host detections

# occupancy models first
# set up detection matrix (host 1, both replicates, both host data)
detections_both_hosts <- data.frame(detections = two_hosts_data$both_hosts)
# set up covariate occupancy frame object (use same covariates)
cov_occu3 <- unmarkedFramePCount(y = detections_both_hosts, siteCovs = site_covs2)
summary(cov_occu3)
# set up various models to compare (just with relative weight gain and fasting period)
pcount_null_3 <- pcount(~ 1 ~1, data = cov_occu3, K = 3)
pcount_m1_3 <- pcount(~ rel_weight_gain ~1, data = cov_occu3, K = 3)
pcount_m2_3 <- pcount(~ rel_weight_gain + fasting_days ~1, data = cov_occu3, K = 3)
# set the fit
fit3 <- fitList('p(.)' = pcount_null_3,
                'p(rel_weight_gain)' = pcount_m1_3,
                'p(rel_weight_gain + fasting_days)' = pcount_m2_3)
modSel(fit3)
# model with fasting days and weight dif not better than null, but still close

# proportional weight predictions
preds3 <- predict(pcount_m2_3, type ="det", newdata = data.frame(rel_weight_gain = seq(0, 1, length = 46),
                                                                fasting_days = seq(0, 45, by =1)))
# plot predictions
ggplot(data = preds3, aes(x = seq(0, 1, length = 46), y = Predicted)) +
  geom_smooth(stat = "smooth") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  theme_classic()

# contour plot as before
xvals2 <- seq(from = min(site_covs2$fasting_days), to = max(site_covs2$fasting_days), length.out = ngrid + 1)
yvals2 <- seq(from = min(site_covs2$rel_weight_gain), to = max(site_covs2$rel_weight_gain), length.out = ngrid + 1)

# generate matrix of occupancy probabilities
x.elev2 <- rep(xvals2, length(xvals2))
y.length2 <- rep(yvals2, each = length(xvals2))

# prediction function
zvals2 <- matrix(predict(pcount_m2_3, type="det", new=data.frame(fasting_days=x.elev2, rel_weight_gain=y.length2))$Predict, nrow=ngrid+1)

# make data suitable for ggplot
contour_data2 <- expand.grid(fasting_days = xvals2, rel_weight_gain = yvals2)
contour_data2$occupancy_prob <- as.vector(zvals2)

# contour plot
ggplot(contour_data2, aes(x = fasting_days, y = rel_weight_gain, z = occupancy_prob)) +
  stat_contour_filled(breaks = seq(0, 1, by = 0.1)) +
  geom_textcontour(linecolour = "#242424", textcolour = "#242424") + 
  labs(x = "Fasting Days", y = "Relative Weight Gain (Host 1)", 
       title = "Both Hosts Predicted Detection Probabilities") +  
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12))











