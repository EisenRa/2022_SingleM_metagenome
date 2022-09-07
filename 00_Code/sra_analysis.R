library(data.table)
library(ggplot2)
library(tidyverse)
library(R.utils)

condensed = fread('zcat ../../../../../../Eisen.Raph/Desktop/sra_20211215.mach3.sdb.condense10pct.mach2.csv.gz')
condensed[1:3]

genome_sizes_dt = data.table(genome_sizes)
genome_sizes_dt[1:3]
nrow(genome_sizes_dt)
genome_sizes_dt[grep('Archaea',rank)]

condensed[, last_taxon := gsub('^.__','',gsub('.*; ', '', taxonomy))]
condensed[, taxonomy := NULL] # delete to save RAM

genome_sizes_dt[,.N,by=rank][N>1][order(-N)][1:4]
genome_sizes_dt[rank=='EX4484-52']

genome_sizes_dt2 = genome_sizes_dt[,.SD[1],by=rank]

estimates = merge(condensed, genome_sizes_dt2, by.x='last_taxon', by.y='rank', all.x=T)[,.(estimated_bp = sum(coverage * genome_size)), by=list(sample)]

estimates[1:3]

metadata = fread('pigz -cd ../../../../../../kingfisher_metadata_some.csv.gz |cut -f1,3,6,8', header=TRUE,sep='\t',quote="")
metadata[1:3]

m = merge(estimates, metadata, by.x='sample', by.y='run')
m[Gbp>2][1:3]

m2 = m[model != 'Illumina Genome Analyzer IIx' & model != 'Illumina Genome Analyzer IIx'][Gbp>2]
m2[, estimated_gbp := estimated_bp/1e9]
m2[, microbial_percent := estimated_gbp/Gbp*100]
m2[1:3]

m2 <- read_delim("00_Code/m2.csv") %>%
  filter(model != "Illumina Genome Analyzer II")

m2 %>%
  group_by(taxon_name) %>%
  summarise(samples = n()) %>%
  arrange(desc(samples), .by_group = TRUE)

types_of_interest = c('human gut metagenome','marine metagenome','soil metagenome', 
                      'food metagenome','human oral metagenome','plant metagenome')

library(glue)

human <- taxa %>% filter(taxon_name == "human gut metagenome") %>% select(n)
marine <- taxa %>% filter(taxon_name == "marine metagenome") %>% select(n)
soil <- taxa %>% filter(taxon_name == "soil metagenome") %>% select(n)
food <- taxa %>% filter(taxon_name == "food metagenome") %>% select(n)
oral <- taxa %>% filter(taxon_name == "human oral metagenome") %>% select(n)
plant <- taxa %>% filter(taxon_name == "plant metagenome") %>% select(n)


names_n <- c(paste0('human gut metagenome', " (n=", {human}, ")"), paste0('marine metagenome', " (n=", {marine}, ")"),
             paste0('soil metagenome', " (n=", {soil}, ")"), paste0('food metagenome', " (n=", {food}, ")"),
             paste0('human oral metagenome', " (n=", {oral}, ")"), paste0('plant metagenome', " (n=", {plant}, ")"))
names(names_n) <- c('human gut metagenome','marine metagenome','soil metagenome', 
                    'food metagenome','human oral metagenome','plant metagenome')


m2 %>%
  filter(taxon_name %in% types_of_interest) %>%
  ggplot(aes(x = microbial_percent)) +
  geom_density(colour = "#e2f0d9", size = 3) +
  facet_wrap(~taxon_name, scales = "free_y", labeller = labeller(taxon_name = names_n)) +
  theme_minimal() +
  xlim(0,100) +
  theme(
    plot.background = element_rect(colour = "transparent", fill = "transparent"),
    axis.text = element_text(size = 14, colour = "white"),
    axis.title = element_text(size = 14, colour = "white"),
    strip.text = element_text(size = 14, face = "bold", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_line(colour = "darkgrey")
  ) +
  xlab("Estimated bacterial/arcaheal percent") +
  ylab("Density")


ggsave("Fig4.png", width=13, height=5)



qplot(data=m2, microbial_percent, geom='density', xlim=c(0,100))

top_taxons = m2[,.N,by=taxon_name][order(-N)][1:20]
top_taxons

theme_set(theme_bw())
types_of_interest = c('human gut metagenome','marine metagenome','soil metagenome','viral metagenome','human oral metagenome','plant metagenome')
qplot(data=m2[taxon_name %in% types_of_interest], microbial_percent, geom='density', xlim=c(0,100), ylab='density', xlab='% of reads that are microbial')+facet_wrap(~taxon_name, scale='free_y')
ggsave('microbial_fraction_estimates.svg', width=11, height=5)

nrow(m2[taxon_name %in% types_of_interest])