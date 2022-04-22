## Code for simulated tests of SingleM's false-positive hit rate on Eukaryotic genomes
export THREADS=40

### Download sequenced Zymo mock communities
#from InStrain paper - cite: https://www.nature.com/articles/s41587-020-00797-0#Sec9

grabseqs sra -t $THREADS SRR12324251 SRR12324252 SRR12324253
mkdir 2_Simulated_reads/Real_zymo_reads
mv 2_Simulated_reads/SRR*.fastq.gz 2_Simulated_reads/Real_zymo_reads

### Spike in simulated human reads into simulated Zymo mock.
#Note, the sequenced Zymo mock samples have:
#SRR12324251 - 8.8 Gbp
#SRR12324252 - 10.8 Gbp
#SRR12324253 - 15.4 Gbp
#I want to randomly subsample these to 0.38 Gbp -- same as the simulated Zymo
for i in 2_Simulated_reads/Real_zymo_reads/*_1.fastq.gz; do
seqtk sample \
-s 1337 $i 1266666 > ${i/_1.fastq.gz/_ss1-2M_1.fastq};
done

for i in 2_Simulated_reads/Real_zymo_reads/*_2.fastq.gz; do
seqtk sample \
-s 1337 $i 1266666 > ${i/_2.fastq.gz/_ss1-2M_2.fastq};
done

pigz -p $THREADS 2_Simulated_reads/Real_zymo_reads/*.fastq

for i in 2_Simulated_reads/Real_zymo_reads/*_ss1-2M_1.fastq.gz; do
cat $i 2_Simulated_reads/Human_reads_3-8M_1.fastq.gz > ${i/_1/_Spiked_1};
done

for i in 2_Simulated_reads/Real_zymo_reads/*_ss1-2M_2.fastq.gz; do
cat $i 2_Simulated_reads/Human_reads_3-8M_2.fastq.gz > ${i/_2/_Spiked_2};
done

#Map reads to Zymo reference genomes
for i in 2_Simulated_reads/Real_zymo_reads/*Spiked_1.fastq.gz; do
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
for i in 2_Simulated_reads/Real_zymo_reads/*Spiked_1.fastq.gz; do
  singlem pipe \
          --singlem_metapackage 0_Database/S3.0.1.metapackage_20211101.smpkg/ \
          --forward $i \
          --reverse ${i/_1.fastq/_2.fastq} \
          --threads $THREADS \
          --include-inserts \
          --otu-table 3_Outputs/$(basename ${i/_1.fastq.gz/.tsv}) \
          --output-extras;
    done
