################################################################################
## Script for estimating bacterial/archaeal fraction of metaegnome using the 
## outputs from SingleM condense and genome_size_from_gtdb.R
## Raphael Eisenhofer July 2022
################################################################################

## Load libraries
library(tidyverse)


  
genome_sizes <- read_delim("3_Outputs/gtdb_mean_genome_sizes.tsv")
singlem_condense <- read_delim("3_Outputs/Simulated_Zymo_Spiked_condense.tsv")

#This regexp will take the last tax rank (detected by '__') from the tax string
str_view(singlem_condense$taxonomy, "[^__]+(?=$)")

#Get lowest taxonomic rank from condense output, match it to the genome size
#table, then create a new column: coverage * genome_size
microbial_fraction_estimate <- singlem_condense %>%
  mutate(lowest_rank = str_extract(taxonomy, "[^__]+(?=$)")) %>%
  inner_join(., genome_sizes, by = c("lowest_rank" = "rank")) %>%
  mutate(estimated_bp = coverage * genome_size,
         estimated_mbp = (coverage * genome_size) / 1000000,
         estimated_gbp = (coverage * genome_size) / 1000000000) 

write_delim(microbial_fraction_estimate, "3_Outputs/Simulated_Zymo_Spiked_microbial_estimate.tsv", delim = "\t")


microbial_fraction_estimate_totals <- microbial_fraction_estimate %>%
  summarise(across(c(estimated_bp, estimated_mbp, estimated_gbp), sum))

write_delim(microbial_fraction_estimate_totals, "3_Outputs/Simulated_Zymo_Spiked_microbial_estimate_totals.tsv", delim = "\t")



