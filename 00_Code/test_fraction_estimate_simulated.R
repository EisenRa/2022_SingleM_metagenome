################################################################################
## Script for estimating bacterial/archaeal fraction of metaegnome using the 
## outputs from SingleM condense and genome_size_from_gtdb.R
## Raphael Eisenhofer July 2022
################################################################################

## load libraries/functions
library(tidyverse)
source("00_Code/estimate_microbial_fraction_function.R")


## Total size of the spiked simulated zymo metagenomes?
sim_zymo_spiked_homo_bp_total <- 757352700 * 2
sim_zymo_spiked_homo_mbp_total <- sim_zymo_spiked_bp_total/1000000
sim_zymo_spiked_arabidopsis_bp_total <- 486523650 * 2
sim_zymo_spiked_arabidopsis_mbp_total <- sim_zymo_spiked_arabidopsis_bp_total/1000000
sim_zymo_spiked_plasmodium_bp_total <- 245668650 * 2
sim_zymo_spiked_plasmodium_mbp_total <- sim_zymo_spiked_plasmodium_bp_total/1000000
  
## Total human bp in spiked simulated zymo metagenome?
sim_zymo_bp <- (187352700 * 2)
sim_zymo_mbp <- sim_zymo_bp/1000000

human_bp <- sim_zymo_spiked_homo_bp_total - sim_zymo_bp
human_mbp <- human_bp/1000000
arabidopsis_bp <- sim_zymo_spiked_arabidopsis_bp_total - sim_zymo_bp
arabidopsis_mbp <- arabidopsis_bp/1000000
plasmodium_bp <- sim_zymo_spiked_plasmodium_bp_total - sim_zymo_bp
plasmodium_mbp <- plasmodium_bp/1000000

## How many mbp of bacterial/arcahaeal DNA are there in our Zymo simulated metagenome?
## This is the ground truth, need to minus the fungi and plasmid seqeunces
ecoli_plasmid <- 110007 * 9
salm_plasmid <- 49572 * 9
staph_plasmid <- 11546 * 15
plasmids <- ecoli_plasmid + salm_plasmid + staph_plasmid
plasmids_mbp <- plasmids/1000000

#N.b. this is # reads in R1 of metagenome * read length * 2
cryptococcus <- 93192 * 150 * 2
saccharomyces <- 42785 * 150 * 2
fungi_bp <- cryptococcus + saccharomyces
fungi_mbp <- fungi_bp/1000000

bacterial_archaeal_ground_truth_bp <- sim_zymo_bp - fungi_bp - plasmids
bacterial_archaeal_ground_truth_mbp <- bacterial_archaeal_ground_truth_bp/1000000

##What's the total number of bact/arch genome copies (coverage) in our simulated zymo metagenomes?
sim_zymo_genome_copy_number <- 99
sim_zymo_spiked_homo_genome_copy_number <- 99
##What about the singlem condense estimates?
sim_zymo_condense <- read_delim("3_Outputs/1_Simulated_zymo/Simulated_Zymo_pipe_condense.tsv")
sim_zymo_spiked_homo_condense <- read_delim("3_Outputs/1_Simulated_zymo/Simulated_Zymo_Spiked_Homo_pipe_condense.tsv")
sim_zymo_spiked_arabidopsis_condense <- read_delim("3_Outputs/1_Simulated_zymo/Simulated_Zymo_Spiked_Arabidopsis_pipe_condense.tsv")
sim_zymo_spiked_plasmodium_condense <- read_delim("3_Outputs/1_Simulated_zymo/Simulated_Zymo_Spiked_Plasmodium_pipe_condense.tsv")
singlem_zymo_genome_copy_number_estimates <- sum(sim_zymo_condense$coverage)
singlem_zymo_spiked_homo_genome_copy_number_estimates <- sum(sim_zymo_spiked_homo_condense$coverage)
singlem_zymo_spiked_arabidopsis_genome_copy_number_estimates <- sum(sim_zymo_spiked_arabidopsis_condense$coverage)
singlem_zymo_spiked_plasmodium_genome_copy_number_estimates <- sum(sim_zymo_spiked_plasmodium_condense$coverage)

#Import the AUGUST improved taxonomic classification approach
sim_zymo_condense_AUG <- read_delim("3_Outputs/1_Simulated_zymo/Simulated_Zymo_AUGUST_pipe_condense.tsv")
zymo_bact_arch_estimate_AUG_bp <- sum(estimate_microbial_fraction(sim_zymo_condense_AUG)$estimated_bp)
zymo_bact_arch_estimate_AUG_mbp <- zymo_bact_arch_estimate_AUG_bp / 1000000

##Different condense parameters
# sim_zymo_condense_10pc_0mtc <- read_delim("3_Outputs/1_Simulated_zymo/Simulated_Zymo_condense_10pc_0mtc.tsv")
# sim_zymo_condense_0pc_0mtc <- read_delim("3_Outputs/1_Simulated_zymo/Simulated_Zymo_condense_0pc_0mtc.tsv")
# sim_zymo_condense_0pc <- read_delim("3_Outputs/1_Simulated_zymo/Simulated_Zymo_condense_0pc.tsv")
# sum(sim_zymo_condense_0pc_0mtc$coverage)
# sum(sim_zymo_condense_0pc$coverage)

##base pair estimates from our method?
zymo_bact_arch_estimate_bp <- sum(estimate_microbial_fraction(sim_zymo_condense)$estimated_bp)
zymo_bact_arch_estimate_mbp <- zymo_bact_arch_estimate_bp / 1000000

zymo_spiked_homo_bact_arch_estimate_bp <- sum(estimate_microbial_fraction(sim_zymo_spiked_homo_condense)$estimated_bp)
zymo_spiked_homo_bact_arch_estimate_mbp <- zymo_spiked_homo_bact_arch_estimate_bp / 1000000

zymo_spiked_arabidopsis_bact_arch_estimate_bp <- sum(estimate_microbial_fraction(sim_zymo_spiked_arabidopsis_condense)$estimated_bp)
zymo_spiked_arabidopsis_bact_arch_estimate_mbp <- zymo_spiked_arabidopsis_bact_arch_estimate_bp / 1000000

zymo_spiked_plasmodium_bact_arch_estimate_bp <- sum(estimate_microbial_fraction(sim_zymo_spiked_plasmodium_condense)$estimated_bp)
zymo_spiked_plasmodium_bact_arch_estimate_mbp <- zymo_spiked_plasmodium_bact_arch_estimate_bp / 1000000


##Pretty plots summarising this:
df <- tibble(
  ground_truth = bacterial_archaeal_ground_truth_mbp,
  our_estimate_species_mean = 313.21,
  our_estimate_genus_mean = 294.8561
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



df2 <- tibble(
  ground_truth = bacterial_archaeal_ground_truth_mbp,
  no_spiking = zymo_bact_arch_estimate_mbp,
  homo_spiked = zymo_spiked_homo_bact_arch_estimate_mbp,
  arabidopsis_spiked = zymo_spiked_arabidopsis_bact_arch_estimate_mbp,
  plasmodium_spiked = zymo_spiked_plasmodium_bact_arch_estimate_mbp,
) %>%
  pivot_longer(everything(), names_to = "method", values_to = "mbp_estimate") %>%
  mutate(percent_estimated = mbp_estimate/332.3029)

write_tsv(df2, "3_Outputs/1_Simulated_zymo/estimate_values.tsv")

df2 %>%
  ggplot(aes(x = method, y = value, fill = method)) +
  geom_bar(stat = "identity") +
  labs(y = "mbp of bacterial/archaeal DNA") +
  theme_classic() +
  theme(
    legend.position = 0
  )


##Why the ~10% discrepancy between ground_truth and estimate?
# Taxonomic specificity of SingleM (i.e. we don't have all species classifications)
# N.B. Lactobacillus fermentum == Limosilactobacillus fermentum on gtdb!!!!
ground_truth_condense <- tibble(
  "__Pseudomonas aeruginosa" = 6,
  "__Escherichia coli" = 9,
  "__Salmonella enterica" = 9,
  "__Limosilactobacillus fermentum" = 21,
  "__Enterococcus faecalis" = 15,
  "__Staphylococcus aureus" = 15,
  "__Listeria monocytogenes" = 14,
  "__Bacillus subtilis" = 10,
) %>%
  pivot_longer(everything(), names_to = "taxonomy", values_to = "coverage")

sum(estimate_microbial_fraction(ground_truth_condense)$estimated_mbp)

bacil_bp <- 4045677
crypto_bp <- 28246122
entero_bp <- 2845392
escher_bp <- 4875441
lacto_bp <- 1905333
lister_bp <- 2992342
pseudo_bp <- 6792330
sacchar_bp <- 12843354
salmo_bp <- 4809318
staph_bp <- 2730326

genome_sizes <- tibble(
  "Bacillus" = 4045677,
  "Cryptococcus" = 28246122,
  "Enterococcus" = 2845392,
  "Escherichia" = 4875441,
  "Lactobacillus" = 1905333,
  "Listeria" = 2992342,
  "Pseudomonas" = 6792330,
  "Saccharomyces" = 12843354,
  "Salmonella" = 4809318,
  "Staphylococcus" = 2730326,
) %>%
  pivot_longer(cols = everything(), names_to = "Genome", values_to = "genome_size")

coverm <- read_delim("3_Outputs/1_Simulated_zymo/Simulated_zymo_Bt2_coverM.tsv")

coverm %>% filter(Genome != "unmapped" & Genome != "Cryptococcus" & Genome != "Saccharomyces") %>%
  select(2) %>%
  summarise(sum = sum(.))

expected_bp <- (bacil_bp * 10)+(entero_bp * 15)+(escher_bp * 9)+(lacto_bp * 21)+(lister_bp * 14)+(pseudo_bp*6)+(salmo_bp * 9)+(staph_bp * 15)

coverm_expected_bp <- coverm %>%
  inner_join(., genome_sizes, by = "Genome") %>%
  mutate(expected_bp_coverm = `Simulated_Zymo_Bt2 Relative Abundance (%)` * genome_size) 

sum(coverm_expected_bp$expected_bp_coverm)

coverm_expected_bp %>%
  filter(Genome != "unmapped" & Genome != "Cryptococcus" & Genome != "Saccharomyces") %>%
  summarise(sum = sum(expected_bp_coverm))


staph <- 136514 * 150 * 2
pseudo <- 135846 * 150 * 2
lacto <- 133373 * 150 * 2
list <- 139642 * 150 * 2
entero <- 142269 * 150 * 2
escher <- 146263 * 150 * 2
salmon <- 144279 * 150 * 2
bacil <- 134855 * 150 * 2

bact <- staph+pseudo+lacto+list+entero+escher+salmon+bacil



