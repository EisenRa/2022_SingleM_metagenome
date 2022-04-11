## Code for simulated tests of SingleM's false-positive hit rate on Eukaryotic genomes

### SingleM database setup
#Download newer SingleM gene package:
wget https://zenodo.org/record/5739612/files/S3.metapackage_20211007.smpkg.tar.gz?download=1
mv S3.metapackage_20211007.smpkg.tar.gz?download=1 S3.metapackage_20211007.smpkg.tar.gz
tar -xvzf S3.metapackage_20211007.smpkg.tar.gz
mv S3.metapackage_20211007.smpkg 0_Database/S3.metapackage_20211007.smpkg
rm S3.metapackage_20211007.smpkg.tar.gz

### Simulating Eukaryotic genome reads
#Download reference genomes (Homo sapien, Arabadopsis thaliana, Plasmodium falciparum):
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/765/GCF_000002765.5_GCA_000002765/GCF_000002765.5_GCA_000002765_genomic.fna.gz
mv *.fna.gz 1_References/

### Create simulated reads from these genomes
#https://github.com/merenlab/reads-for-assembly
#N.B., takes uncompressed fastas as input
pigz -d 1_References/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gen-paired-end-reads 00_Code/Human_config.ini
pigz -p 40 1_References/GCF_000001405.39_GRCh38.p13_genomic.fna

pigz -d 1_References/GCF_000001735.4_TAIR10.1_genomic.fna.gz
gen-paired-end-reads 00_Code/Arabadopsis_config.ini
pigz -p 40 1_References/GCF_000001735.4_TAIR10.1_genomic.fna

pigz -d 1_References/GCF_000002765.5_GCA_000002765_genomic.fna.gz
gen-paired-end-reads 00_Code/Plasmodium_config.ini
pigz -p 40 1_References/GCF_000002765.5_GCA_000002765_genomic.fna

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

singlem pipe \
--singlem_metapackage S3.metapackage_20211007.smpkg/ \
--forward 2_Simulated_reads/Arabadopsis_reads-R1.fastq.gz \
--reverse 2_Simulated_reads/Arabadopsis_reads-R2.fastq.gz \
--threads 40 \
--include-inserts \
--otu-table Arabidopsis_singleM_inserts \
--output-extras

singlem pipe \
--singlem_metapackage S3.metapackage_20211007.smpkg/ \
--forward 2_Simulated_reads/Human_reads-R1.fastq.gz \
--reverse 2_Simulated_reads/Human_reads-R2.fastq.gz \
--threads 40 \
--include-inserts \
--otu-table Human_singleM_inserts \
--output-extras
