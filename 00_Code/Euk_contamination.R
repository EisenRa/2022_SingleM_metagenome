################################################################################
### Script for testing/plotting eukaryotic false positives in SingleM
###
################################################################################

library(tidyverse)
library(scales)

arabidopsis_pipe <- read_delim("3_Outputs/0_Host_genome_false_positive/Arabidopsis_singleM_inserts.csv")
plasmodium_pipe <- read_delim("3_Outputs/0_Host_genome_false_positive/Plasmodium_singleM.csv")
homo_pipe <- read_delim("3_Outputs/0_Host_genome_false_positive/Human_singleM.csv")

arabidopsis_condense <- read_delim("3_Outputs/0_Host_genome_false_positive/Arabidopsis_singleM_inserts_condense.csv")
plasmodium_condense <- read_delim("3_Outputs/0_Host_genome_false_positive/Plasmodium_singleM_condense.csv")
homo_condense <- read_delim("3_Outputs/0_Host_genome_false_positive/Human_singleM_condense.csv")

arabidopsis_n_hits <- sum(arabidopsis_pipe$num_hits)
plasmodium_n_hits <- sum(plasmodium_pipe$num_hits)
homo_n_hits <- sum(homo_pipe$num_hits)

arabidopsis_condense_sum <- sum(arabidopsis_condense$coverage)
plasmodium_condense_sum <- sum(plasmodium_condense$coverage)
homo_condense_sum <- sum(homo_condense$coverage)

#Euk genome size information
arabidopsis_genome_size_mb <- 119.669
plasmodium_genome_size_mb <- 23.3328
homo_genome_size_mb <- 3298.43

arabidopsis_mb_total <- arabidopsis_genome_size_mb * 5
plasmodium_mb_total <- plasmodium_genome_size_mb * 5
homo_mb_total <- homo_genome_size_mb * 2

#N.B. these are paired reads
arabidopsis_n_sim_reads <- 1994473
plasmodium_n_sim_reads <- 388773
homo_n_sim_reads <- 21813601

false_positives <- tibble(
  "arabidopsis_cov" = arabidopsis_condense_sum,
  "plasmodium_cov" = plasmodium_condense_sum,
  "homo_cov" = homo_condense_sum,
  "arabidopsis_mb_total" = arabidopsis_mb_total,
  "plasmodium_mb_total" = plasmodium_mb_total,
  "homo_mb_total" = homo_mb_total
) %>%
  pivot_longer(cols = everything(),
               names_to = c("genome", ".value"),
               names_sep = "_",
               values_to = c("cov", "mb")) %>%
  mutate(cov_normalised = cov / mb)
  
  

#Plotto
false_positives %>%
  ggplot(aes(x = genome, y = cov_normalised, fill = genome)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(
    legend.position = "none"
  ) +
  ylab("false positive coverage of bacteria/archaea / Mbp host genome")

ggsave("3_Outputs/Euk_false_positive_percent_barplot.png")
