## Code for creating simulated Zymo mock communities

### Download Zymo mock communities reference genomes
wget https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.refseq.v2.zip


### Create simulated reads from these genomes
#https://github.com/merenlab/reads-for-assembly
#N.B., takes uncompressed fastas as input
pigz -d 1_References/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gen-paired-end-reads 00_Code/Simulated_Zymo.ini
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
