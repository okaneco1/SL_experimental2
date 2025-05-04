# Analysis of Host Switching Data for Experiment 2

options(scipen = 999)

# libraries
library(tidyverse)
library(readxl)
library(patchwork)
library(RColorBrewer)
library(gridExtra)
library(pscl)
library(ggpubr)
library(MuMIn)
library(MASS)


# import data
host_data <- read_csv("host_switch_data_averaged.csv")
cm_data_avg <- read_excel("submission1_averaged_community_matrix.xlsx")
cm_data_reps <- read_csv("host_switch_data_replicates.csv")

# create data frame for total sample reads (average)
seq_count <- cm_data_avg %>%
  filter(grepl("T23", sample)) %>%
  # specific samples to remove
  filter(!grepl("T23_30", sample) & !grepl("T23_18", sample)) %>%
  dplyr::select(sample, Salmonidae_unclassified_mean, Salvelinus_namaycush_mean, 
         Catostomus_commersonii_mean, Petromyzontidae_unclassified_mean,
         Salvelinus_unclassified_mean, Catostomidae_unclassified_mean, 
         Catostomus_unclassified_mean) %>% # selecting lake trout, white sucker, and sea lamprey 
  mutate(total_reads = rowSums(across(-1)),
         all_trout = rowSums(across(c(Salmonidae_unclassified_mean, Salvelinus_namaycush_mean)))) %>%
  rename(white_sucker = Catostomus_commersonii_mean)

#prepare for bar chart visualization
seq_count_long <- pivot_longer(seq_count, 
                               cols = c(white_sucker, all_trout),
                               names_to = "species",
                               values_to = "read_count") %>%
  dplyr::select(sample, species, read_count)
seq_count_long$species <- factor(seq_count_long$species, levels = c("all_trout", "white_sucker"))

# adding the fasting days for x axis labels
seq_metadata <- host_data %>%
  dplyr::select(sample, fasting_days) %>%
  left_join(
    cm_data_reps %>% 
      filter(!is.na(sample)) %>%
      dplyr::select(sample, weight_gain_2),
    by = "sample"
  )


seq_count_long_days <- seq_count_long %>%
  left_join(seq_metadata, by = "sample") %>%
  mutate(fasting_days = na_if(fasting_days, "n/a")) %>%
  filter(!is.na(fasting_days)) %>%
  filter(weight_gain_2 > 0) %>%
  arrange(fasting_days)

# set labels for fasting days
fasting_labels <- seq_count_long_days %>%
  distinct(sample, fasting_days) %>%
  deframe() 

# refactor fasting days for proper order
seq_count_long_days <- seq_count_long_days %>%
  mutate(fasting_days = as.numeric(fasting_days)) %>%
  mutate(sample = factor(sample, levels = unique(sample[order(fasting_days)])))

  
# plotting
ggplot(seq_count_long_days, aes(x = sample, y = read_count, fill = species)) +
  geom_col() +
  scale_fill_manual(
    values = c("all_trout" = "#4c8a75", "white_sucker" = "#ad9be0"),
    labels = c("white_sucker" = "White Sucker", "all_trout" = "Lake Trout"),
    name = "Species") +
  #scale_x_discrete(labels = fasting_labels) +
  labs(x = "Fasting Days",
       y = "Host Fish Sequence Reads") +
  theme_minimal(base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_line(color = alpha("gray", 0.1)),
    panel.grid.minor = element_line(color = alpha("gray", 0.1))
  )




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
    dplyr::select(fasting_days, host1_det_rep1, host1_det_rep2) %>%
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

write.csv(host_data_reps_clean, file = "host_data_reps_clean.csv")




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

# adding a host 2 read count column
host_data_reps_clean$host2_read_count_rep1 <- NA
host_data_reps_clean$host2_read_count_rep2 <- NA

for (i in 1:nrow(host_data_reps_clean)) {
  if (host_data_reps_clean$host_2[i] == "Lake Trout") {
    host_data_reps_clean$host2_read_count_rep1[i] <- host_data_reps_clean$lake_trout_reads_rep1[i]
    host_data_reps_clean$host2_read_count_rep2[i] <- host_data_reps_clean$lake_trout_reads_rep2[i]
  }
  else if (host_data_reps_clean$host_2[i] == "White Sucker") {
    host_data_reps_clean$host2_read_count_rep1[i] <- host_data_reps_clean$white_sucker_reads_rep1[i]
    host_data_reps_clean$host2_read_count_rep2[i] <- host_data_reps_clean$white_sucker_reads_rep2[i]
  }
}

### add average host2 read count column

# viewing data
hist(log(host_data_reps_clean$host1_read_count_rep1))
hist(host_data_reps_clean$host1_read_count_rep2)
# reshape for plot
long_data <- host_data_reps_clean %>%
  pivot_longer(cols = c(host1_read_count_rep1, host1_read_count_rep2),
               names_to = "replicate",
               values_to = "read_count")

# plot each variable against read count (for each rep)
ggplot(long_data, aes(x = rel_weight_gain, y = read_count, color = replicate)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()
ggplot(long_data, aes(x = fasting_days, y = read_count, color = replicate)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, span = 1) + # span increased to desensitize fit
  geom_jitter() +
  theme_minimal()
ggplot(long_data, aes(x = host1_species, y = read_count, color = replicate)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()

# adding "fasting since host 1" table
host_data_reps_clean <- host_data_reps_clean %>%
  mutate(fasting_days_1 = fasting_days + 7)

# Model Comparisons

# will first combine rep1 and rep2 read counts so that they are both considered
host_data_reps_clean <- mutate(host_data_reps_clean, host1_read_count_avg = rowMeans(cbind(host1_read_count_rep1,
                                                                                           host1_read_count_rep2)))
host_data_reps_clean <- mutate(host_data_reps_clean, host2_read_count_avg = rowMeans(cbind(host2_read_count_rep1,
                                                                                           host2_read_count_rep2)))

# setting up models (need to round for models that accept count data)
lm1 <- lm((host1_read_count_avg) ~ rel_weight_gain + host1_species + fasting_days * fasting_days_1 + days_attached_1, data = host_data_reps_clean)
negbin_lm <- glm.nb(round(host1_read_count_avg) ~ rel_weight_gain + host1_species + fasting_days * fasting_days_1 + days_attached_1, na.action = na.fail, data = host_data_reps_clean)
zinb_lm <- zeroinfl(round(host1_read_count_avg) ~ rel_weight_gain + days_attached_1 + fasting_days + host1_species | 1, 
                    data = host_data_reps_clean, 
                    dist = "negbin")
# AIC comparisoons
AIC(lm1)
AIC(negbin_lm)
AIC(zinb_lm)

summary(negbin_lm)

# negative binomial shows highest AIC again, so can use that

# use dredge to check all combinations
dredge_table <- dredge(negbin_lm)
dredge_df <- data.frame(dredge_table)
write.csv(dredge_df, file = "negative_binomial_predictor_table_exp2.csv", row.names = FALSE)

# null model is lowest AIC, but days attached is close
1-exp(-0.2847) # 0.2477601



summary(negbin_lm)

#test
filtered_days_data <- host_data_reps_clean %>%
  filter(days_attached_1 > 4)

# quick visual of averaged read count (host 1) vs fasting days
ggplot(host_data_reps_clean, aes(x = days_attached_1, y = host1_read_count_avg)) +
  geom_jitter(width = 0.1, height=0.3) +
  geom_smooth(method = "lm", se = FALSE, span = 1) + # span increased to desensitize fit
  theme_minimal()



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



# Comparing Read Counts of Host 1 and Host 2

# averaged between two replicates
# reshape data
read_count_host_long <- host_data_reps_clean %>%
  filter(weight_gain_2 > 0) %>% # remove fish that did not feed on second host
  dplyr::select(host1_read_count_avg, host2_read_count_avg) %>%
  pivot_longer(cols = c(host1_read_count_avg, host2_read_count_avg),
               names_to = "host",
               values_to = "read_count_avg")
  
# rename host column variables
read_count_host_long$host <- recode(read_count_host_long$host,
                                    "host1_read_count_avg" = "host 1",
                                    "host2_read_count_avg" = "host 2")

# anova test
host_reads_lm <- lm(read_count_avg ~ host, data = read_count_host_long)
anova_result <- anova(host_reads_lm)

# make plot
ggplot(read_count_host_long, aes(x = host, y = read_count_avg)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3) +
  labs(
       x = "Host Number",
       y = "Read Count Average") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12), 
    axis.title.x = element_text(face = "bold"),          
    axis.title.y = element_text(face = "bold"))

summary_stats <- read_count_host_long %>%
  group_by(host) %>%
  summarize(
    mean_read_count = mean(read_count_avg),
    median_read_count = median(read_count_avg),
    min_read_count = min(read_count_avg),
    max_read_count = max(read_count_avg),
    sd_read_count = sd(read_count_avg)
  )


# ------------------------------------------------------------------
# looking at removing if hosts did NOT feed on a second host

# this is because if there is no second host feeding, detection of 
# an feeding history is nullified, as only a single host was fed on
# ------------------------------------------------------------------

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
  labs(x = "Fasting Days", y = "Detection Rate (Host 1)", 
       title = "A)",
       color = "Replicate") +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("#87CEEB", "#228B22"),
                     labels = c("Replicate 1", "Replicate 2")) +
  geom_hline(yintercept = 0, color = "#5c5e5c", linewidth = 0.5) +   
  geom_vline(xintercept = 0, color = "#5c5e5c", linewidth = 0.5) +
  scale_x_continuous(breaks = c(0, 5, 10, 20, 30)) +  
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        legend.position = "inside",
        legend.title = element_text(face = "bold"),
        legend.position.inside = c(0.2, 0.7),
        legend.background = element_rect(fill = alpha('white', 0.5), color = NA))

# NOTE: warnings do not seem to make a difference compared to non-jittered graph


#-----------------------------------------------------
# Descriptive Statistics
# Two Host Samples
#-----------------------------------------------------
# set up detection data frame, along with column for at least one det between reps
detection_df <- two_hosts_data %>%
  dplyr::select(sample, host1_det_rep1, host1_det_rep2, fasting_days) %>%
  mutate(host1_det_both = ifelse(host1_det_rep1 == 1 | host1_det_rep2 == 1, 1, 0))

detection_df_stats <- detection_df %>%
  summarize(
    detection_percentage = mean(host1_det_both),
    detection_sd = sd(host1_det_both)
    )

# graph to summarize data
fasting_stats <- data.frame(
  fasting_days = unique(two_hosts_data$fasting_days),
  detection_total = NA,
  samples = NA,
  det_proportion = NA
) %>%
  arrange(fasting_days)

# fill in data
for (i in 1:nrow(fasting_stats)){
  # set up detection data frame
  rep_df <- two_hosts_data %>%
    dplyr::select(fasting_days, host1_det_rep1, host1_det_rep2) %>%
    filter(fasting_days == fasting_stats$fasting_days[i]) %>%
    mutate(any_host_det = ifelse(host1_det_rep1 | host1_det_rep2 == 1, 1, 0))
  # add in detection total
  fasting_stats$detection_total[i] <- sum(rep_df$any_host_det)
  # add in sample number
  fasting_stats$samples[i] <- nrow(rep_df)
  # add in detection proportion
  fasting_stats$det_proportion[i] <- round(fasting_stats$detection_total[i]/fasting_stats$samples[i], 3)
}

fasting_stats
  


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
  labs(x = "Fasting Days", y = "Detection Rate (Host 2)", 
       title = "B)",
       color = "Replicate") +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("#87CEEB", "#228B22"),
                     labels = c("Replicate 1", "Replicate 2")) +
  geom_hline(yintercept = 0, color = "#5c5e5c", linewidth = 0.5) +   
  geom_vline(xintercept = 0, color = "#5c5e5c", linewidth = 0.5) +
  scale_x_continuous(breaks = c(0, 5, 10, 20, 30)) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
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
both_hosts_fasting_plot <- ggplot(two_hosts_data, aes(x = fasting_days, y = both_hosts)) +
  geom_jitter(width = 0.5, height = 0.02, color = "#003582") +
  labs(x = "Fasting Days", y = "Detection Rate (Both Hosts)", 
       title = "C)") +
  geom_smooth(method = "loess", se = FALSE, color = "#003582") +
  geom_hline(yintercept = 0, color = "#5c5e5c", linewidth = 0.5) +   
  geom_vline(xintercept = 0, color = "#5c5e5c", linewidth = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10))

# plot all together
grid.arrange(host1_two_hosts_fasting_plot, 
             host2_two_hosts_fasting_plot,
             both_hosts_fasting_plot, ncol=1)





