################################################################################
## Testing out the estimate microbial fraction function on CAMI2 datasets 
## 
## Raphael Eisenhofer July 2022
################################################################################

## load libraries/functions
library(tidyverse)
source("00_Code/estimate_microbial_fraction_function.R")


## Metagenome stats
marine0_paired_reads <- 16647395
marine0_bp <- 2497109250 * 2
marine0_mbp <- marine0_bp / 1000000
marine0_gbp <- marine0_mbp / 1000
marine0_genome_copies_new <- 1427.699658

plant0_paired_reads <- 16627054
plant0_bp <- 2494058100 * 2
plant0_mbp <- plant0_bp / 1000000
plant0_gbp <- plant0_mbp / 1000
plant0_genome_copies_new <- 658.226053

smadness0_paired_reads <- 6655023
smadness0_bp <- 998253450 * 2
smadness0_mbp <- smadness0_bp / 1000000
smadness0_gbp <- smadness0_mbp / 1000
smadness0_genome_copies_new <- 512.4978345

## Import condense outputs
marine0_condense <- read_delim("3_Outputs/3_CAMI2/marine_0_condense.tsv")
plant0_condense <- read_delim("3_Outputs/3_CAMI2/plant_0_condense.tsv")
smadness0_condense <- read_delim("3_Outputs/3_CAMI2/strain_madness_0_condense.tsv")

# SingleM estimated genome copy numbers?
marine0_singlem_estimated_coverage <- sum(marine0_condense$coverage)
marine0_singlem_performance <- marine0_singlem_estimated_coverage/marine0_genome_copies_new

plant0_singlem_estimated_coverage <- sum(plant0_condense$coverage)
plant0_singlem_performance <- plant0_singlem_estimated_coverage/plant0_genome_copies_new

smadness0_singlem_estimated_coverage <- sum(smadness0_condense$coverage)
smadness0_singlem_performance <- smadness0_singlem_estimated_coverage/smadness0_genome_copies_new

## Get bacterial/archaeal fraction estimates
marine0_estimate_bp <- sum(estimate_microbial_fraction(marine0_condense)$estimated_bp)
marine0_estimate_bp / 1000000000

plant0_estimate_bp <- sum(estimate_microbial_fraction(plant0_condense)$estimated_bp)
plant0_estimate_bp / 1000000000

smadness0_estimate_bp <- sum(estimate_microbial_fraction(smadness0_condense)$estimated_bp)
smadness0_estimate_bp / 1000000000


gt <- read_delim("3_Outputs/3_CAMI2/coverage_new0.tsv", col_names = c("genome", "coverage"))

gt_plasmids <- gt %>% 
  filter(str_detect(genome, "^RNODE")) %>%
  separate(., genome, sep = "_", into = c("1", "2", "3", "length", "5", "6"), convert = TRUE) %>%
  mutate(bp = length * coverage)

gt_bact_arch <- gt %>%
  filter(!str_detect(genome, "^RNODE"))

sum(gt_bact_arch$coverage)
