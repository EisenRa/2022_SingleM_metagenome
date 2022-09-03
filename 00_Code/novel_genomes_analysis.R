

library(tidyverse)
library(janitor)

genome_stats <- read_delim("1_References/Novel_genomes/quality_report.tsv") %>%
  mutate(Name = str_replace(Name, ".fna", ""))

df <- read_delim("1_References/Novel_genomes/gtdbtk.bac120.summary.tsv", na = "N/A") %>%
  left_join(., genome_stats, by = c("user_genome" = "Name"))

df_filtered_placement_ani <- df %>%
  filter(Completeness > 80 & Contamination < 10 & closest_placement_ani < 95) %>%
  mutate(user_genome = str_replace(user_genome, "_genomic", ""))

df_filtered_no_reference <- df %>%
  filter(Completeness > 80 & Contamination < 10 & is.na(fastani_reference)) %>%
  mutate(user_genome = str_replace(user_genome, "_genomic", ""))

write_tsv(df_filtered_placement_ani %>% 
            select(user_genome) %>% 
            mutate(user_genome = word(user_genome, 1,2, sep = "_")),
          "1_References/Novel_genomes/genomes_placement_ani.tsv",
          col_names = FALSE)

write_tsv(df_filtered_no_reference %>% 
            select(user_genome) %>% 
            mutate(user_genome = word(user_genome, 1,2, sep = "_")),
          "1_References/Novel_genomes/genomes_no_reference.tsv",
          col_names = FALSE)

