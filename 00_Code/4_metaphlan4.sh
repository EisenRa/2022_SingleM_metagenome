#Issue with quality headers -- read_fastx errors out because the "+" has text after it.
for i in ../2022_SingleM_metagenome/2_Simulated_reads/Simulated_Zymo*R2.fastq.gz;
do \
zcat $i | sed 's/+.*/+/g' > $(basename ${i/.gz/});
done

pigz -p 40 *.fastq

for i in *_R1.fastq.gz;
do \
metaphlan "$i",${i/R1/R2} \
--bowtie2db ../../0_DBs/ \
--bowtie2out ${i/.fastq.gz/.bz2} \
--nproc 40 \
--input_type fastq \
--unclassified_estimation \
-o ${i/.fastq.gz/.txt};
done
