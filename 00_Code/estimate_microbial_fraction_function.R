################################################################################
## Script for estimating bacterial/archaeal fraction of metaegnome using the 
## outputs from SingleM condense and genome_size_from_gtdb.R
## Raphael Eisenhofer July 2022
################################################################################

## Load libraries
library(tidyverse)

estimate_microbial_fraction <- function(condense_file) {

  
genome_sizes <- read_delim("3_Outputs/gtdb_r207_mean_genome_sizes.tsv")

#Get lowest taxonomic rank from condense output, match it to the genome size
#table, then create a new column: coverage * genome_size
microbial_fraction_estimate <- condense_file %>%
  mutate(lowest_rank = word(taxonomy, -1, sep = "; ")) %>%
  inner_join(., genome_sizes, by = c("lowest_rank" = "rank")) %>%
  mutate(estimated_bp = coverage * genome_size,
         estimated_mbp = (coverage * genome_size) / 1000000,
         estimated_gbp = (coverage * genome_size) / 1000000000) 

return(microbial_fraction_estimate)

}

