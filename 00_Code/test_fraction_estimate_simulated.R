################################################################################
## Script for estimating bacterial/archaeal fraction of metaegnome using the 
## outputs from SingleM condense and genome_size_from_gtdb.R
## Raphael Eisenhofer July 2022
################################################################################

## load libraries/functions
library(tidyverse)
source("00_Code/estimate_microbial_fraction_function.R")


## Total size of human spiked simulated zymo metagenome?
sim_zymo_spiked_bp_total <- 757352700 * 2
sim_zymo_spiked_mbp_total <- sim_zymo_spiked_bp_total/1000000

## Total human bp in spiked simulated zymo metagenome?
sim_zymo_bp <- (187352700 * 2)
sim_zymo_mbp <- sim_zymo_bp/1000000

human_bp <- sim_zymo_spiked_bp_total - sim_zymo_bp
human_mbp <- human_bp/1000000

## How many mbp of bacterial/arcahaeal DNA are there in our Zymo simulated metagenome?
## This is the ground truth, need to minus the fungi and plasmid seqeunces
ecoli_plasmid <- 110007 * 9
salm_plasmid <- 49572 * 9
staph_plasmid <- 11546 * 15
plasmids <- ecoli_plasmid + salm_plasmid + staph_plasmid

#N.b. this is # reads in R1 of metagenome * read length * 2
cryptococcus <- 93192 * 150 * 2
saccharomyces <- 42785 * 150 * 2
fungi_bp <- cryptococcus + saccharomyces

bacterial_archaeal_ground_truth_bp <- sim_zymo_bp - fungi_bp - plasmids
bacterial_archaeal_ground_truth_mbp <- bacterial_archaeal_ground_truth_bp/1000000

##What's the total number of bact/arch genome copies (coverage) in our simulated zymo metagenomes?
sim_zymo_genome_copy_number <- 99
sim_zymo_spiked_genome_copy_number <- 99
##What about the singlem condense estimates?
sim_zymo_condense <- read_delim("3_Outputs/1_Simulated_zymo/Simulated_Zymo_condense_10pc.tsv")
sim_zymo_human_spiked_condense <- read_delim("3_Outputs/1_Simulated_zymo/Simulated_Zymo_Spiked_condense.tsv")
singlem_zymo_genome_copy_number_estimates <- sum(sim_zymo_condense$coverage)
singlem_zymo_spiked_genome_copy_number_estimates <- sum(sim_zymo_spiked_condense$coverage)

##Different condense parameters
sim_zymo_condense_10pc_0mtc <- read_delim("3_Outputs/1_Simulated_zymo/Simulated_Zymo_condense_10pc_0mtc.tsv")
sim_zymo_condense_0pc_0mtc <- read_delim("3_Outputs/1_Simulated_zymo/Simulated_Zymo_condense_0pc_0mtc.tsv")
sim_zymo_condense_0pc <- read_delim("3_Outputs/1_Simulated_zymo/Simulated_Zymo_condense_0pc.tsv")
sum(sim_zymo_condense_0pc_0mtc$coverage)
sum(sim_zymo_condense_0pc$coverage)

sum(estimate_microbial_fraction(sim_zymo_condense_0pc_0mtc)$estimated_mbp)
sum(estimate_microbial_fraction(sim_zymo_condense_0pc)$estimated_mbp)


##base pair estimates from our method?
zymo_bact_arch_estimate_bp <- sum(estimate_microbial_fraction(sim_zymo_condense)$estimated_bp)
zymo_bact_arch_estimate_mbp <- zymo_bact_arch_estimate_bp / 1000000

zymo_spiked_bact_arch_estimate_bp <- sum(estimate_microbial_fraction(sim_zymo_spiked_condense)$estimated_bp)
zymo_spiked_bact_arch_estimate_mbp <- bact_arch_estimate_bp / 1000000


##Pretty plots summarising this:
df <- tibble(
  ground_truth = bacterial_archaeal_ground_truth_mbp,
  our_estimate = zymo_spiked_bact_arch_estimate_mbp
) %>%
  pivot_longer(everything(), names_to = "method")

df %>%
  ggplot(aes(x = method, y = value, fill = method)) +
  geom_bar(stat = "identity") +
  labs(y = "mbp of bacterial/archaeal DNA") +
  theme_classic() +
  theme(
    legend.position = 0
  )

ggsave("3_Outputs/bacterial_archaeal_estimate_barplot.png")


##Why the ~10% discrepancy between ground_truth and estimate?
# Taxonomic specificity of SingleM (i.e. we don't have all species classifications)
ground_truth_condense <- tibble(
  "__Pseudomonas aeruginosa" = 6,
  "__Escherichia coli" = 9,
  "__Salmonella enterica" = 9,
  "__Lactobacillus fermentum" = 21,
  "__Enterococcus faecalis" = 15,
  "__Staphylococcus aureus" = 15,
  "__Listeria monocytogenes" = 14,
  "__Bacillus subtilis" = 10,
) %>%
  pivot_longer(everything(), names_to = "taxonomy", values_to = "coverage")

sum(estimate_microbial_fraction(ground_truth_condense)$estimated_mbp)


staph <- 136514 * 150 * 2
pseudo <- 135846 * 150 * 2
lacto <- 133373 * 150 * 2
list <- 139642 * 150 * 2
entero <- 142269 * 150 * 2
escher <- 146263 * 150 * 2
salmon <- 144279 * 150 * 2
bacil <- 134855 * 150 * 2

bact <- staph+pseudo+lacto+list+entero+escher+salmon+bacil



