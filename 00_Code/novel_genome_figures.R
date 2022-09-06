################################################################################
## Script for estimating bacterial/archaeal fraction of metaegnome using the 
## outputs from SingleM condense and genome_size_from_gtdb.R
## Raphael Eisenhofer July 2022
################################################################################

## load libraries/functions
library(tidyverse)
source("00_Code/estimate_microbial_fraction_function.R")

## Total size of the simulated novel metagenomes:
novel_ani_bp <- 515953050 * 2
novel_ani_mbp <- novel_ani_bp / 1000000
novel_ref_bp <- 1057906950 * 2
novel_ref_mbp <- novel_ref_bp / 1000000

##Import singlem condense results
novel_ani_condense <- read_delim("3_Outputs/5_Novel_genomes/Novel_genomes_placement_ani_AUGUST_pipe_condense.tsv")
novel_ref_condense <- read_delim("3_Outputs/5_Novel_genomes/Novel_genomes_no_reference_AUGUST_pipe_condense.tsv")

novel_ani_estimate_bp <- sum(estimate_microbial_fraction(novel_ani_condense)$estimated_bp)
novel_ani_estimate_mbp <- novel_ani_estimate_bp / 1000000
novel_ani_estimate_percent <- novel_ani_estimate_mbp / novel_ani_mbp

novel_ref_estimate_bp <- sum(estimate_microbial_fraction(novel_ref_condense)$estimated_bp)
novel_ref_estimate_mbp <- novel_ref_estimate_bp / 1000000
novel_ref_estimate_percent <- novel_ref_estimate_mbp / novel_ref_mbp

##Tibblize!
df <- tibble(
  `novel community 1` = novel_ani_estimate_percent,
  `novel community 2` = novel_ref_estimate_percent
) %>%
  pivot_longer(everything(), names_to = "method", values_to = "percent_estimated")

## Poster plot!

figure3_col <- c("novel community 1" = "#00b0f0", "novel community 2" = "#e2f0d9")

df %>%
  ggplot(aes(x = method, 
             y = percent_estimated * 100, 
             fill = method)
  ) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 100, colour = "yellow", linetype="dashed", size = 2) +
  labs(y = "Percent of bacterial DNA estimated") +
  theme_classic() +
  scale_fill_manual(values = figure3_col) +
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

ggsave("PosterFigure3.png", width = 5, height = 6, units = "in", bg = 'transparent')



