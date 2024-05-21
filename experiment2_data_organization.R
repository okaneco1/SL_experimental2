# Merging Metabarcoding and Sample Data: Experiment 2 - Host-Swtiching Trials


# libraries
library(tidyverse)
library(readxl)

# import data
cm_data_avg <- read_excel("submission1_averaged_community_matrix.xlsx")
cm_data_rep1 <- read_excel("submission1A_community_matrix.xls", range = "A1:AW506")
cm_data_rep2 <- read_excel("submission1B_community_matrix.xls")


#---------- Data Organization

# community matrix data for averages and separate replications will be 
# processed individually to allow for analysis of both

# -------
# function for processing
process_cm_data <- function(data, lake_trout_cols, white_sucker_cols) {
  data %>%
    filter(grepl("T23", sample)) %>%
    filter(!grepl("T23_30", sample)) %>%
    mutate(
      lake_trout_reads = rowSums(select(., all_of(lake_trout_cols))),
      white_sucker_reads = rowSums(select(., all_of(white_sucker_cols)))
    ) %>%
    select(sample, lake_trout_reads, white_sucker_reads)
}

# --------
# process data for averaged community matrix
cm_data_avg <- process_cm_data(cm_data_avg,
                lake_trout_cols = c("Salmonidae_unclassified_mean",
                                    "Salvelinus_namaycush_mean",
                                    "Salvelinus_unclassified_mean"),
                white_sucker_cols = c("Catostomus_commersonii_mean",
                                      "Catostomidae_unclassified_mean",
                                      "Catostomus_unclassified_mean"))

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
                                                     "Catostomus_unclassified"))
cm_data_rep2 <- process_cm_data(cm_data_rep2,
                                lake_trout_cols = c("Salmonidae_unclassified",
                                                    "Salvelinus_namaycush",
                                                    "Salvelinus_unclassified"),
                                white_sucker_cols = c("Catostomus_commersonii",
                                                      "Catostomidae_unclassified",
                                                      "Catostomus_unclassified"))

# making names the same
cm_data_rep1$sample <- sub("\\.12S", "", cm_data_rep1$sample)
cm_data_rep2$sample <- sub("_R\\.12S", "", cm_data_rep2$sample)
# reorder
cm_data_rep1 <- cm_data_rep1[order(cm_data_rep1[[1]], decreasing = FALSE), ]
cm_data_rep2 <- cm_data_rep2[order(cm_data_rep2[[1]], decreasing = FALSE), ]




