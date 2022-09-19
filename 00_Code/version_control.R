################################################################################
## Escaping version control hell 
## 
## Raphael Eisenhofer Sept 2022
################################################################################

## load libraries/functions
library(tidyverse)
source("00_Code/estimate_microbial_fraction_function.R")


## Metagenome stats
marine0_paired_reads <- 16647395
marine0_groundtruth_bp <- 2497105500 * 2
marine0_groundtruth_plasmids_bp <- 56371500 * 2
marine0_groundtruth_bacteria_bp <- marine0_groundtruth_bp - marine0_groundtruth_plasmids_bp
marine0_groundtruth_mbp <- marine0_groundtruth_bacteria_bp / 1000000
marine0_groundtruth_gbp <- marine0_groundtruth_mbp / 1000
marine0_genome_copies_new <- 1427.699658

plant0_paired_reads <- 16627054
plant0_groundtruth_bp <- 2494058100 * 2
plant0_groundtruth_mbp <- plant0_groundtruth_bp / 1000000
plant0_groundtruth_gbp <- plant0_groundtruth_mbp / 1000
plant0_genome_copies_new <- 658.226053

smadness0_paired_reads <- 6655023
smadness0_groundtruth_bp <- 998253450 * 2
smadness0_groundtruth_mbp <- smadness0_groundtruth_bp / 1000000
smadness0_groundtruth_gbp <- smadness0_groundtruth_mbp / 1000
smadness0_genome_copies_new <- 512.4978345

## Import condense outputs
marine0_CAMISIM_condense <- read_delim("3_Outputs/3_CAMI2/marine_0_condense.tsv")
plant0_CAMISIM_condense <- read_delim("3_Outputs/3_CAMI2/plant_0_condense.tsv")
smadness0_CAMISIM_condense <- read_delim("3_Outputs/3_CAMI2/strain_madness_0_condense.tsv")

#August results:
marine0_CAMISIM_august_condense <- read_delim("3_Outputs/3_CAMI2/marine_0_AUGUST_pipe_condense.tsv")
plant0_CAMISIM_august_condense <- read_delim("3_Outputs/3_CAMI2/plant_0_AUGUST_pipe_condense.tsv")
smadness0_CAMISIM_august_condense <- read_delim("3_Outputs/3_CAMI2/strain_madness_0_AUGUST_pipe_condense.tsv")

#August results (home-made simulated metagenomes)
marine0_HM_august_condense <- read_delim("3_Outputs/3_CAMI2/CAMI2_marine0_HomeMade_AUGUST_pipe_condense.tsv")
plant0_HM_august_condense <- read_delim("3_Outputs/3_CAMI2/CAMI2_plant0_HomeMade_AUGUST_pipe_condense.tsv")
smadness0_HM_august_condense <- read_delim("3_Outputs/3_CAMI2/CAMI2_strain0_HomeMade_AUGUST_pipe_condense.tsv")

#Sept deployment test
marine0_HM_september_condense <- read_delim("3_Outputs/1_Simulated_zymo/test_sept_marine_condense.tsv")
marine0_HM_september_estimate_bp <- sum(estimate_microbial_fraction(marine0_HM_september_condense)$estimated_bp)
marine0_HM_september_estimate_mbp <- marine0_HM_september_estimate_bp/1000000


# SingleM estimated genome copy numbers?
marine0_CAMISIM_singlem_estimated_coverage <- sum(marine0_CAMISIM_condense$coverage)
marine0_CAMISIM_singlem_performance <- marine0_CAMISIM_singlem_estimated_coverage/marine0_genome_copies_new

plant0_CAMISIM_singlem_estimated_coverage <- sum(plant0_CAMISIM_condense$coverage)
plant0_CAMISIM_singlem_performance <- plant0_CAMISIM_singlem_estimated_coverage/plant0_genome_copies_new

smadness0_CAMISIM_singlem_estimated_coverage <- sum(smadness0_CAMISIM_condense$coverage)
smadness0_CAMISIM_singlem_performance <- smadness0_CAMISIM_singlem_estimated_coverage/smadness0_genome_copies_new

#August results
marine0_CAMISIM_august_singlem_estimated_coverage <- sum(marine0_CAMISIM_august_condense$coverage)
marine0_CAMISIM_august_singlem_performance <- marine0_CAMISIM_august_singlem_estimated_coverage/marine0_genome_copies_new

plant0_CAMISIM_august_singlem_estimated_coverage <- sum(plant0_CAMISIM_august_condense$coverage)
plant0_CAMISIM_august_singlem_performance <- plant0_CAMISIM_august_singlem_estimated_coverage/plant0_genome_copies_new

smadness0_CAMISIM_august_singlem_estimated_coverage <- sum(smadness0_CAMISIM_august_condense$coverage)
smadness0_CAMISIM_august_singlem_performance <- smadness0_CAMISIM_august_singlem_estimated_coverage/smadness0_genome_copies_new

#August results (home-made simulated metagenomes)
marine0_HM_august_singlem_estimated_coverage <- sum(marine0_HM_august_condense$coverage)
marine0_HM_august_singlem_performance <- marine0_HM_august_singlem_estimated_coverage/marine0_genome_copies_new

plant0_HM_august_singlem_estimated_coverage <- sum(plant0_HM_august_condense$coverage)
plant0_HM_august_singlem_performance <- plant0_HM_august_singlem_estimated_coverage/plant0_genome_copies_new

smadness0_HM_august_singlem_estimated_coverage <- sum(smadness0_HM_august_condense$coverage)
smadness0_HM_august_singlem_performance <- smadness0_HM_august_singlem_estimated_coverage/smadness0_genome_copies_new


## Get bacterial/archaeal fraction estimates
marine0_CAMISIM_estimate_bp <- sum(estimate_microbial_fraction(marine0_CAMISIM_condense)$estimated_bp)
marine0_CAMISIM_estimate_mbp <- marine0_CAMISIM_estimate_bp/1000000

plant0_CAMISIM_estimate_bp <- sum(estimate_microbial_fraction(plant0_CAMISIM_condense)$estimated_bp)
plant0_CAMISIM_estimate_mbp <- plant0_CAMISIM_estimate_bp/1000000

smadness0_CAMISIM_estimate_bp <- sum(estimate_microbial_fraction(smadness0_CAMISIM_condense)$estimated_bp)
smadness0_CAMISIM_estimate_mbp <- smadness0_CAMISIM_estimate_bp/1000000

# August results
marine0_CAMISIM_august_estimate_bp <- sum(estimate_microbial_fraction(marine0_CAMISIM_CAMISIM_august_condense)$estimated_bp)
marine0_CAMISIM_august_estimate_mbp <- marine0_CAMISIM_august_estimate_bp/1000000

plant0_CAMISIM_august_estimate_bp <- sum(estimate_microbial_fraction(plant0_CAMISIM_august_condense)$estimated_bp)
plant0_CAMISIM_august_estimate_mbp <- plant0_CAMISIM_august_estimate_bp/1000000

smadness0_CAMISIM_august_estimate_bp <- sum(estimate_microbial_fraction(smadness0_CAMISIM_august_condense)$estimated_bp)
smadness0_CAMISIM_august_estimate_mbp <- smadness0_CAMISIM_august_estimate_bp/1000000

#August results (home-made simulated metagenomes)
marine0_HM_august_estimate_bp <- sum(estimate_microbial_fraction(marine0_HM_august_condense)$estimated_bp)
marine0_HM_august_estimate_mbp <- marine0_HM_august_estimate_bp/1000000

plant0_HM_august_estimate_bp <- sum(estimate_microbial_fraction(plant0_HM_august_condense)$estimated_bp)
plant0_HM_august_estimate_mbp <- plant0_HM_august_estimate_bp/1000000

smadness0_HM_august_estimate_bp <- sum(estimate_microbial_fraction(smadness0_HM_august_condense)$estimated_bp)
smadness0_HM_august_estimate_mbp <- smadness0_HM_august_estimate_bp/1000000


## Plots

df_all <- tibble(
  marine0_HM_august_estimate_mbp, smadness0_HM_august_estimate_mbp,
  marine0_CAMISIM_august_estimate_mbp, smadness0_CAMISIM_august_estimate_mbp,
  marine0_groundtruth_mbp, smadness0_groundtruth_mbp
) %>%
  pivot_longer(everything(), names_to = "environment", values_to = "estimate") %>%
  mutate(origin = case_when(str_detect(environment, "HM") ~ "home_made",
                            str_detect(environment, "CAMISIM") ~ "CAMISIM",
                            str_detect(environment, "groundtruth") ~ "ground_truth"),
         community = case_when(str_detect(environment, "marine") ~ "marine",
                               str_detect(environment, "smad") ~ "strain-madness")
  )

df_all %>%
  ggplot(aes(x = reorder(environment, estimate), 
             y = estimate, 
             fill = factor(origin, levels = c("CAMISIM", "home_made", "ground_truth")))) +
  geom_bar(stat = "identity") +
  facet_wrap(~community, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  ) +
  ylab("mbp")

ggsave("3_Outputs/3_CAMI2/SingleM_CAMI2_HomeBrew.pdf", width = 10, height = 10, unit = "in")



df <- tibble(
  marine0_HM_august_estimate_mbp, smadness0_HM_august_estimate_mbp,
  marine0_groundtruth_mbp, smadness0_groundtruth_mbp
) %>%
  pivot_longer(everything(), names_to = "environment", values_to = "estimate") %>%
  mutate(source = case_when(str_detect(environment, "HM") ~ "home_made",
                            str_detect(environment, "CAMISIM") ~ "CAMISIM",
                            str_detect(environment, "groundtruth") ~ "ground_truth"),
         community = case_when(str_detect(environment, "marine") ~ "marine",
                               str_detect(environment, "smad") ~ "strain-madness")
  )

df %>%
  ggplot(aes(x = environment, y = estimate, fill = source)) +
  geom_bar(stat = "identity") +
  facet_wrap(~community, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  ylab("mbp")

ggsave("3_Outputs/3_CAMI2/SingleM_CAMI2.pdf", width = 10, height = 10, unit = "in")


#PERCENTAGES

df_all_perc <- tibble(
  marine0_HM_august_singlem_performance, smadness0_HM_august_singlem_performance,
  marine0_CAMISIM_august_singlem_performance, smadness0_CAMISIM_august_singlem_performance
) %>%
  pivot_longer(everything(), names_to = "environment", values_to = "estimate") %>%
  mutate(origin = case_when(str_detect(environment, "HM") ~ "home_made",
                            str_detect(environment, "CAMISIM") ~ "CAMISIM"),
         community = case_when(str_detect(environment, "marine") ~ "marine",
                               str_detect(environment, "smad") ~ "strain-madness")
  )

df_all_perc %>%
  ggplot(aes(x = reorder(environment, estimate), 
             y = estimate, 
             fill = factor(origin, levels = c("CAMISIM", "home_made")))) +
  geom_bar(stat = "identity") +
  facet_wrap(~community, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  ) +
  ylab("mbp")

ggsave("3_Outputs/3_CAMI2/SingleM_CAMI2_HomeBrew_PERCENTAGE.pdf", width = 10, height = 10, unit = "in")



## Poster figure

df4 <- tibble(marine0_HM_august_singlem_performance, 
              smadness0_HM_august_singlem_performance) %>%
  pivot_longer(everything(), names_to = "environment", values_to = "estimate") %>%
  mutate(origin = case_when(str_detect(environment, "HM") ~ "home_made",
                            str_detect(environment, "CAMISIM") ~ "CAMISIM"),
         community = case_when(str_detect(environment, "marine") ~ "marine",
                               str_detect(environment, "smad") ~ "strain-madness")
  )

figure2_col <- c("marine" = "#002060", "strain-madness" = "RED")


df4 %>%
  ggplot(aes(x = community, 
             y = estimate * 100, 
             fill = community)
  )+
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 100, colour = "yellow", linetype="dashed", size = 2) +
  labs(y = "Percent of bacterial DNA estimated") +
  theme_classic() +
  scale_fill_manual(values = figure2_col) +
  coord_cartesian(expand = FALSE) +
  theme(
    legend.position = 0,
    axis.text = element_text(size = 15, colour = "white"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18, colour = "white"),
    plot.background = element_rect(fill = "transparent", colour = 'transparent'),
    panel.background = element_rect(fill = 'transparent', colour = 'transparent'),
    axis.ticks = element_line(colour = "#e2f0d9"),
    axis.line = element_line(colour = "#e2f0d9")
  ) 

ggsave("PosterFigure2.png", width = 5, height = 6, units = "in", bg = 'transparent')











#Check coverage of plasmids

gt_plasmids <- gt %>% 
  filter(str_detect(genome, "^RNODE")) %>%
  separate(., genome, sep = "_", into = c("1", "2", "3", "length", "5", "6"), convert = TRUE) %>%
  mutate(bp = length * coverage)

gt_bact_arch <- gt %>%
  filter(!str_detect(genome, "^RNODE"))

sum(gt_bact_arch$coverage)
