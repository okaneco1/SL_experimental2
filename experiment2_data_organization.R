# Merging Metabarcoding and Sample Data: Experiment 2 - Host-Swtiching Trials


# libraries
library(tidyverse)
library(readxl)

# import data
cm_data_avg <- read_excel("submission1_averaged_community_matrix.xlsx")
cm_data_rep1 <- read_excel("submission1A_community_matrix.xls", range = "A1:AW506")
cm_data_rep2 <- read_excel("submission1B_community_matrix.xls")
dissection_data <- read_excel("dissection_data_2023_updated2.xlsx", range = "A1:M70")
host1_data <- read_csv("obj3_cleaned_data_host1_updated.csv")
host2_data <- read_csv("obj3_cleaned_data_host2_updated.csv")



#---------- Data Organization

# community matrix data for averages and separate replications will be 
# processed individually to allow for analysis of both

# -------
# function for processing
process_cm_data <- function(data, lake_trout_cols, white_sucker_cols, sea_lamprey_cols) {
  data %>%
    filter(grepl("T23", sample)) %>%
    filter(!grepl("T23_30", sample)) %>%
    mutate(
      lake_trout_reads = rowSums(select(., all_of(lake_trout_cols))), 
      white_sucker_reads = rowSums(select(., all_of(white_sucker_cols))),
      sea_lamprey_reads = rowSums(select(., all_of(sea_lamprey_cols)))
    ) %>%
    select(sample, lake_trout_reads, white_sucker_reads, sea_lamprey_reads)
}

# --------
# process data for averaged community matrix
cm_data_avg <- process_cm_data(cm_data_avg,
                lake_trout_cols = c("Salmonidae_unclassified_mean",
                                    "Salvelinus_namaycush_mean",
                                    "Salvelinus_unclassified_mean"),
                white_sucker_cols = c("Catostomus_commersonii_mean",
                                      "Catostomidae_unclassified_mean",
                                      "Catostomus_unclassified_mean"),
                sea_lamprey_cols = c("Petromyzontidae_unclassified_mean",
                                     "Ichthyomyzon_fossor_mean",
                                     "Lampetra_appendix_mean"))

# -------
# ensure "sample" is first column for replicate datasets
colnames(cm_data_rep1)[1] <- "sample"
colnames(cm_data_rep2)[1] <- "sample"

# process data for individual replication community matrix data
cm_data_rep1 <- process_cm_data(cm_data_rep1,
                               lake_trout_cols = c("Salmonidae_unclassified",
                                                   "Salvelinus_namaycush",
                                                   "Salvelinus_unclassified"),
                               white_sucker_cols = c("Catostomus_commersonii",
                                                     "Catostomidae_unclassified",
                                                     "Catostomus_unclassified"),
                               sea_lamprey_cols = c("Petromyzontidae_unclassified",
                                                    "Ichthyomyzon_fossor",
                                                    "Lampetra_appendix"))
cm_data_rep2 <- process_cm_data(cm_data_rep2,
                                lake_trout_cols = c("Salmonidae_unclassified",
                                                    "Salvelinus_namaycush",
                                                    "Salvelinus_unclassified"),
                                white_sucker_cols = c("Catostomus_commersonii",
                                                      "Catostomidae_unclassified",
                                                      "Catostomus_unclassified"),
                                sea_lamprey_cols = c("Petromyzontidae_unclassified",
                                                     "Ichthyomyzon_fossor",
                                                     "Lampetra_appendix"))

# making names the same
cm_data_rep1$sample <- sub("\\.12S", "", cm_data_rep1$sample)
cm_data_rep2$sample <- sub("_R\\.12S", "", cm_data_rep2$sample)
# reorder
cm_data_rep1 <- cm_data_rep1[order(cm_data_rep1[[1]], decreasing = FALSE), ]
cm_data_rep2 <- cm_data_rep2[order(cm_data_rep2[[1]], decreasing = FALSE), ]


# ------- combine host 1 and host 2 data
host1_columns <- data.frame(tag = host1_data$`Lamprey Tag (color)`,
                            host = host1_data$`Host Species`,
                            weight_gain = host1_data$`Change in Weight(g)`,
                            days_attached = host1_data$`Days Attached`)
host2_columns <- data.frame(tag = host2_data$`Lamprey Tag (color)`,
                            host = host2_data$`Host species`,
                            weight_gain = host2_data$`Change in Weight(g)`,
                            days_attached = host2_data$`Days Attached`)

# combine both
host_data <- full_join(host1_columns, host2_columns, by = "tag",
                       suffix = c("_1", "_2"))

#--- checking for duplicates

# may be some duplicate tag numbers, so let's look into those
duplicate_tags_host1 <- host1_columns %>%
  group_by(tag) %>%
  filter(n() > 1) %>%
  pull(tag) %>%
  unique()

duplicate_tags_host2 <- host2_columns %>%
  group_by(tag) %>%
  filter(n() > 1) %>%
  pull(tag) %>%
  unique()

duplicate_tags_host1
duplicate_tags_host2

# previous run throughs had duplicates, but updates to the data now leave
# duplicates for only host 2 with two "No Tag" entries (no actual IDs)
# so, good to go here



#------- combine tag ID data with lamprey feeding/host data
# change column name in host data
colnames(host_data)[1] <- "lamprey"

# start by linking lamprey tags to Tube IDs
tag_tube_data <- data.frame(lamprey = dissection_data$`Lamprey ID`,
                            tube = dissection_data$`Tube ID`,
                            fasting_days = dissection_data$`Fast (days)`)
tag_tube_data <- tag_tube_data[-30, ] # this was a 2022 sample

# join tube IDs with lamprey tag numbers
full_host_data <- host_data %>%
  full_join(tag_tube_data, by = "lamprey")

# seem to be some tubes that are missing corresponding lamprey samples
# can filter these out to list them:
full_host_data %>%
  filter(is.na(host_1) & is.na(host_2)) %>%
  pull(tube) 
# "T2023_13": died early, did not attach
# "T2023_27": may need to discard
# "T2023_60": may need to discard
# "T2023_65": died early, did not attach

# in previous runs there were duplicate tubes as well
duplicate_tubes <- full_host_data %>%
  group_by(lamprey) %>%
  filter(n() > 1) %>%
  pull(lamprey) %>%
  unique()
duplicate_tubes

full_host_data %>%
  filter(lamprey %in% duplicate_tubes)
# now, only "no tag" tubes that were detected previously as well
# will just discard these from the analysis


#--------- Aligning Host Data with Sequence Data
# first re-label tubes to match with the host data
full_host_data$tube <- sub("T20(\\d{2}_\\d+)", "T\\1", full_host_data$tube)


# can combine with REPLICATE SPECIFIC DATA

# start with combining replicate read data
cm_both_reps <- full_join(cm_data_rep1, cm_data_rep2, by = "sample", suffix = c("_rep1", "_rep2"))
# match up sample name regex patterns (ensure that number has two digits)
cm_both_reps$sample <- sub("T23_(\\d)$", "T23_0\\1", cm_both_reps$sample)
# join
full_reps_data <- full_host_data %>%
  rename(sample = tube) %>%
  full_join(cm_both_reps, by = "sample")
# write out
write_csv(full_reps_data, "host_switch_data_replicates.csv")











#-------- old code
# first making a data frame linking host order to lamprey
host_order <- data.frame(lamprey = rep(NA, nrow(host_data)), # empty data frame
                         host_1 = NA, host_2 = NA)
for (i in 1:nrow(host_data)) {
  if (host_data$host_number[i] == 1) {
    host_order$lamprey[i] <- host_data$`Lamprey Tag (color)`[i]
    host_order$host_1[i] <- host_data$`Host Species`[i]
  }
  if (host_data$host_number[i] == 2) {
    host_order$lamprey[i] <- host_data$`Lamprey Tag (color)`[i]
    host_order$host_2[i] <- host_data$`Host Species`[i]
  }
}
# group hosts with lamprey
host_order <- host_order %>% 
  group_by(lamprey) %>%
  summarize(
    host_1 = first(na.omit(host_1)),
    host_2 = first(na.omit(host_2))
  ) %>%
  ungroup()




# can align host order data with sequencing data (AVERAGED COMMUNITY MATRIX)
# full join
full_data <- full_host_data %>%
  rename(sample = tube) %>%
  full_join(cm_data_avg, by = "sample")
# reorder
full_data <- full_data[,c(4,2:3,6:7,5)]
# arrange sample order
full_data <- full_data %>% arrange(sample)

# write out as own data file
write_csv(full_data, "host_switch_data_averaged.csv")
