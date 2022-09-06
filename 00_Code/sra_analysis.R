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

qplot(data=m2, microbial_percent, geom='density', xlim=c(0,100))

top_taxons = m2[,.N,by=taxon_name][order(-N)][1:20]
top_taxons

theme_set(theme_bw())
types_of_interest = c('human gut metagenome','marine metagenome','soil metagenome','viral metagenome','human oral metagenome','plant metagenome')
qplot(data=m2[taxon_name %in% types_of_interest], microbial_percent, geom='density', xlim=c(0,100), ylab='density', xlab='% of reads that are microbial')+facet_wrap(~taxon_name, scale='free_y')
ggsave('microbial_fraction_estimates.svg', width=11, height=5)

nrow(m2[taxon_name %in% types_of_interest])