## Code for creating simulated Zymo mock communities
export THREADS=40

### Download Zymo mock communities reference genomes
wget https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.refseq.v2.zip
unzip ZymoBIOMICS.STD.refseq.v2.zip

### Create simulated reads from these genomes
#https://github.com/merenlab/reads-for-assembly
#N.B., takes uncompressed fastas as input
pigz -d 1_References/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gen-paired-end-reads 00_Code/Simulated_Zymo.ini
pigz -p $THREADS 1_References/GCF_000001405.39_GRCh38.p13_genomic.fna

pigz -p $THREADS 2_Simulated_reads/*.fastq

### Map simulated reads to reference genomes using Bowtie2
cat 1_References/ZymoBIOMICS.STD.refseq.v2/Genomes/*.fasta > 1_References/ZymoCatted.fna
pigz -p $THREADS 1_References/ZymoCatted.fna

bowtie2-build \
--threads $THREADS \
1_References/ZymoCatted.fna.gz 1_References/ZymoCatted.fna.gz

bowtie2 \
--threads $THREADS \
-x 1_References/ZymoCatted.fna.gz \
-1  \
-2  \
--seed 1337 \
| samtools sort -@ $THREADS -o Simulated_Zymo_Bt2.bam -

coverm genome \
        -b 2_Simulated_reads/*.bam \
        -s _ \
        -m relative_abundance \
        -t $THREADS \
        --min-covered-fraction 0 \
        > 3_Outputs/Simulated_zymo_Bt2_coverM.tsv

### Run SingleM pipe on the simulated reads
singlem pipe \
--singlem_metapackage S3.metapackage_20211007.smpkg/ \
--forward 2_Simulated_reads/Plasmodium_reads-R1.fastq.gz \
--reverse 2_Simulated_reads/Plasmodium_reads-R2.fastq.gz \
--threads $THREADS \
--include-inserts \
--otu-table Plasmodium_singleM_inserts \
--output-extras
