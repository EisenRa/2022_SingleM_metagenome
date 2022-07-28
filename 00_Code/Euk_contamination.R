################################################################################
### Script for testing/plotting eukaryotic false positives in SingleM
###
################################################################################

library(tidyverse)
library(scales)

arabidopsis_pipe <- read_delim("3_Outputs/Arabidopsis_singleM_inserts.csv")
plasmodium_pipe <- read_delim("3_Outputs/Plasmodium_singleM.csv")
homo_pipe <- read_delim("3_Outputs/Human_singleM.csv")

arabidopsis_condense <- read_delim("3_Outputs/Arabidopsis_singleM_inserts_condense.csv")
plasmodium_condense <- read_delim("3_Outputs/Plasmodium_singleM_condense.csv")
homo_condense <- read_delim("3_Outputs/Human_singleM_condense.csv")

arabidopsis_n_hits <- nrow(arabidopsis_pipe)
plasmodium_n_hits <- nrow(plasmodium_pipe)
homo_n_hits <- nrow(homo_pipe)

#Euk genome size information
arabidopsis_genome_size_mb <- 119.669
plasmodium_genome_size_mb <- 23.3328
homo_genome_size_mb <- 3298.43

arabidopsis_mb_total <- arabidopsis_genome_size_mb * 5
plasmodium_mb_total <- plasmodium_genome_size_mb * 5
homo_mb_total <- homo_genome_size_mb * 2

arabidopsis_n_sim_reads <- 1994473 * 2
plasmodium_n_sim_reads <- 388773 * 2
homo_n_sim_reads <- 21813601 * 2

n_hits_df <- tibble(arabidopsis_n_hits, plasmodium_n_hits, homo_n_hits) %>%
  pivot_longer(., everything(), names_to = "genome", values_to = "n_hits", names_pattern = "(.*)_n_hits") 

n_reads_df <- tibble(arabidopsis_n_sim_reads, plasmodium_n_sim_reads, homo_n_sim_reads) %>%
  pivot_longer(., everything(), names_to = "genome", values_to = "n_reads", names_pattern = "(.*)_n_sim_reads")

false_positives <- n_hits_df %>%
  left_join(., n_reads_df, by = ("genome"))

false_positives <- false_positives %>%
  mutate(percent = n_hits / n_reads)

#Plotto
false_positives %>%
  ggplot(aes(x = genome, y = percent, fill = genome)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = percent) +
  theme_classic() +
  theme(
    legend.position = "none"
  )
