###
###
###


## Load libraries
library(tidyverse)


## Simulated Zymo mocks
##Import
sim_zymo <- read_delim("3_Outputs/Simulated_Zymo.tsv")
sim_zymo_spiked <- read_delim("3_Outputs/Simulated_Zymo_Spiked.tsv")
sim_zymo_coverm <- read_delim("3_Outputs/Simulated_zymo_Bt2_coverM.tsv")


## Simulated Zymo
#Remove Eukaryota hits, remove outliers
sim_zymo_by_gene <- sim_zymo %>%
  filter(str_detect(taxonomy, ".Eukaryota.", negate = TRUE)) %>%
  group_by(gene) %>%
  summarise(n_hits = sum(num_hits), total_cov = sum(coverage)) %>%
  filter(total_cov >= 70 & total_cov <= 130)

#Get average of coverage across all filtered markers
sim_zymo_by_gene %>%
  summarise(mean = mean(total_cov))

#Expected coverage of simulated Zymo mock = 100 -- i.e. 100 genome copies

#Coverage from mapping sim reads to ref genomes (Bowtie2)
sim_zymo_coverm %>%
  filter(Genome != "unmapped") %>%
  summarise(sum = sum(`Simulated_Zymo_Bt2 Relative Abundance (%)`))


## Spiked (with sim human reads) simulated Zymo mock
#Remove Eukaryota hits, remove outliers
sim_zymo_by_gene <- sim_zymo_spiked %>%
  filter(str_detect(taxonomy, ".Eukaryota.", negate = TRUE)) %>%
  group_by(gene) %>%
  summarise(n_hits = sum(num_hits), total_cov = sum(coverage)) %>%
  filter(total_cov >= 70 & total_cov <= 130)


#Coverage from mapping sim reads to ref genomes (Bowtie2)
sim_zymo_coverm %>%
  filter(Genome != "unmapped") %>%
  summarise(sum = sum(`Simulated_Zymo_Spiked_Bt2 Relative Abundance (%)`))



## Real world Zymo mocks (spiked with sim human reads)
real_zymo1 <- read_delim("3_Outputs/SRR12324251_ss1-2M_Spiked.tsv")
real_zymo2 <- read_delim("3_Outputs/SRR12324252_ss1-2M_Spiked.tsv")
real_zymo3 <- read_delim("3_Outputs/SRR12324253_ss1-2M_Spiked.tsv")
real_zymo_coverm <- read_delim("3_Outputs/Real_zymo_Bt2_coverM.tsv")


#Remove Eukaryota hits, remove outliers
real_zymo1_by_gene <- real_zymo1 %>%
  filter(str_detect(taxonomy, ".Eukaryota.", negate = TRUE)) %>%
  group_by(gene) %>%
  summarise(n_hits = sum(num_hits), total_cov = sum(coverage)) %>%
  filter(total_cov >= 70 & total_cov <= 130)

#Get average of coverage across all filtered markers
real_zymo1_by_gene %>%
  summarise(mean = mean(total_cov))

#Coverage from mapping sim reads to ref genomes (Bowtie2)
real_zymo_coverm %>%
  filter(Genome != "unmapped") %>%
  summarise(sum = sum(`SRR12324251_ss1-2M_Spiked_Bt2 Relative Abundance (%)`))
