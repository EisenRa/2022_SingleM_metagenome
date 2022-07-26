##
##
##

## Load libraries
library(tidyverse)
library(janitor)


## Download data from GTDB and decompress
download.file("https://data.gtdb.ecogenomic.org/releases/release207/207.0/ar53_metadata_r207.tar.gz",
              "0_Database/ar53_metadata_r207.tar.gz")
download.file("https://data.gtdb.ecogenomic.org/releases/release207/207.0/bac120_metadata_r207.tar.gz",
              "0_Database/bac120_metadata_r207.tar.gz")
untar("0_Database/ar53_metadata_r207.tar.gz", exdir = "0_Database/")
untar("0_Database/bac120_metadata_r207.tar.gz", exdir = "0_Database/")

## Import data
gtdb_bac_metadata <- read_tsv("0_Database/bac120_metadata_r207.tsv")
gtdb_arc_metadata <- read_tsv("0_Database/ar53_metadata_r207.tsv")

## Bind bac/arc metadata, filter by minimnum/max completeness/contamination
## Also split out taxonomy strink into columns by taxonomic rank
gtdb_combined_metadata <- rbind(gtdb_arc_metadata, gtdb_bac_metadata) %>%
  filter(checkm_completeness > 80 & checkm_contamination < 10) %>%
  mutate(gtdb_taxonomy = str_replace_all(gtdb_taxonomy, ".__", "")) %>%
  separate(., col = gtdb_taxonomy,
           sep = ";",
           into = c("domain", "phylum", "class", "order", "family", "genus", "species")
           )

################################################################################
## Basic stats + plots

gtdb_combined_metadata %>%
  summarise(mean = mean(genome_size), median = median(genome_size),
            max = max(genome_size), min = min(genome_size),
            n_phyla = n_distinct(phylum), n_family = n_distinct(family),
            n_species = n_distinct(species))

phyla_g_sizes <- gtdb_combined_metadata %>%
  group_by(phylum) %>%
  summarise(mean = mean(genome_size), std = sd(genome_size),
            min = min(genome_size), max = max(genome_size))

family_g_sizes <- gtdb_combined_metadata %>%
  group_by(family) %>%
  summarise(mean = mean(genome_size), std = sd(genome_size),
            min = min(genome_size), max = max(genome_size))

genus_g_sizes <- gtdb_combined_metadata %>%
  group_by(genus) %>%
  summarise(mean = mean(genome_size), std = sd(genome_size),
            min = min(genome_size), max = max(genome_size))


gtdb_combined_metadata %>%
  ggplot(., aes(x = domain, y = genome_size)) +
  geom_boxplot()


################################################################################
## Mean genome sizes per taxonomic rank
genome_size_domain <- gtdb_combined_metadata %>%
  group_by(domain) %>%
  summarise(mean = mean(genome_size)) %>%
  mutate(rank = "domain") %>%
  rename("name" = domain)

genome_size_phylum <- gtdb_combined_metadata %>%
  group_by(phylum) %>%
  summarise(mean = mean(genome_size)) %>%
  mutate(rank = "phylum") %>%
  rename("name" = phylum)

genome_size_class <- gtdb_combined_metadata %>%
  group_by(class) %>%
  summarise(mean = mean(genome_size)) %>%
  mutate(rank = "class") %>%
  rename("name" = class)

genome_size_order <- gtdb_combined_metadata %>%
  group_by(order) %>%
  summarise(mean = mean(genome_size)) %>%
  mutate(rank = "order") %>%
  rename("name" = order)

genome_size_family <- gtdb_combined_metadata %>%
  group_by(family) %>%
  summarise(mean = mean(genome_size)) %>%
  mutate(rank = "family") %>%
  rename("name" = family)

genome_size_genus <- gtdb_combined_metadata %>%
  group_by(genus) %>%
  summarise(mean = mean(genome_size)) %>%
  mutate(rank = "genus") %>%
  rename("name" = genus)

genome_size_species <- gtdb_combined_metadata %>%
  group_by(species) %>%
  summarise(mean = mean(genome_size)) %>%
  mutate(rank = "species") %>%
  rename("name" = species)

full_genome_size_table <- rbind(genome_size_domain, genome_size_phylum,
                                genome_size_class, genome_size_order,
                                genome_size_family, genome_size_genus,
                                genome_size_species)

write_tsv(full_genome_size_table, "gtdb_mean_genome_sizes.tsv")
