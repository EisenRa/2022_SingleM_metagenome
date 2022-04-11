## Code for simulated tests of SingleM's false-positive hit rate on Eukaryotic genomes
export THREADS=40

### Download sequenced Zymo mock communities
#from InStrain paper - cite: https://www.nature.com/articles/s41587-020-00797-0#Sec9

grabseqs sra -t $THREADS SRR12324251 SRR12324252 SRR12324253
mv *.fastq.gz 2_Simulated_reads/Real_zymo_reads

### Create simulated reads from these genomes
#https://github.com/merenlab/reads-for-assembly
#N.B., takes uncompressed fastas as input
pigz -d 1_References/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gen-paired-end-reads 00_Code/Human_config_spikein.ini
pigz -p $THREADS 1_References/GCF_000001405.39_GRCh38.p13_genomic.fna

pigz -p $THREADS 2_Simulated_reads/*.fastq


#Map reads to Zymo reference genomes
for i in 2_Simulated_reads/Real_zymo_reads/*_1.fastq.gz; do
  bowtie2 \
          --threads $THREADS \
          -x 1_References/ZymoCatted.fna.gz \
          -1 $i \
          -2 ${i/_1.fastq/_2.fastq} \
          --seed 1337 \
          | samtools sort -@ $THREADS -o ${i/_1.fastq.gz/_Bt2.bam} -;
    done

coverm genome \
        -b 2_Simulated_reads/Real_zymo_reads/*.bam \
        -s _ \
        -m relative_abundance \
        -t $THREADS \
        --min-covered-fraction 0 \
        > 3_Outputs/Real_zymo_Bt2_coverM.tsv


### Run SingleM pipe on the Zymo reads
for i in 2_Simulated_reads/Real_zymo_reads/*_1.fastq.gz; do
  singlem pipe \
          --singlem_metapackage S3.metapackage_20211007.smpkg/ \
          --forward $i \
          --reverse ${i/_1.fastq/_2.fastq} \
          --threads $THREADS \
          --include-inserts \
          --otu-table 3_Outputs/$(basename ${i/_1.fastq.gz/}) \
          --output-extras;
    done
