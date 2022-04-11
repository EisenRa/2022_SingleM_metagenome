## Code for simulated tests of SingleM's false-positive hit rate on Eukaryotic genomes

### Download sequenced Zymo mock communities
#from InStrain paper - cite: https://www.nature.com/articles/s41587-020-00797-0#Sec9

SRR12324251
SRR12324252
SRR12324253


### Create simulated reads from these genomes
#https://github.com/merenlab/reads-for-assembly
#N.B., takes uncompressed fastas as input
pigz -d 1_References/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gen-paired-end-reads 00_Code/Human_config_spikein.ini
pigz -p 40 1_References/GCF_000001405.39_GRCh38.p13_genomic.fna

pigz -p 40 2_Simulated_reads/*.fastq



### Run SingleM pipe on the simulated reads
singlem pipe \
--singlem_metapackage S3.metapackage_20211007.smpkg/ \
--forward 2_Simulated_reads/Plasmodium_reads-R1.fastq.gz \
--reverse 2_Simulated_reads/Plasmodium_reads-R2.fastq.gz \
--threads 40 \
--include-inserts \
--otu-table Plasmodium_singleM_inserts \
--output-extras
