
cd 1_Reference/Novel_genomes

#Pull accessions from NCBI using ncbi-genome-download:
ncbi-genome-download -s genbank bacteria,archaea -F fasta -p 40 -A genomes_placement_ani.tsv
mkdir novel_ani && mv genbank/bacteria/*/*.gz novel_ani/
ncbi-genome-download -s genbank bacteria,archaea -F fasta -p 40 -A genomes_no_reference.tsv
mkdir novel_ref && mv genbank/bacteria/*/*.gz novel_ref/

# Some of these assemblies have fragments shorter than 1,000 bp, so remove them
module load tools
module load jdk/18.0.1 bbmap/38.90

for i in novel_ani/*.fna.gz;
  do reformat.sh \
  in=$i \
  out=${i/.gz/} \
  minlength=1000;
done

for i in novel_ref/*.fna.gz;
  do reformat.sh \
  in=$i \
  out=${i/.gz/} \
  minlength=1000;
done

## Create coverage profiles to feed into gen-paired-reads
# Add general header information for gen-paired-reads
echo -e '[general]\noutput_sample_name = 2_Simulated_reads/Novel_genomes_placement_ani\ninsert_size = 30\ninsert_size_std = 1\nshort_read_length = 150\nerror_rate = 0.005\n' > Novel_ani.txt
echo -e '[general]\noutput_sample_name = 2_Simulated_reads/Novel_genomes_no_reference\ninsert_size = 30\ninsert_size_std = 1\nshort_read_length = 150\nerror_rate = 0.005\n' > Novel_ref.txt

# Create coverage profiles
for i in novel_ani/*.fna;
  do echo "[$PWD/$i]" && echo "coverage = 3";
done >> novel_ani_coverage.txt

for i in novel_ref/*.fna;
  do echo "[$PWD/$i]" && echo "coverage = 3";
done >> novel_ref_coverage.txt

cat Novel_ani.txt novel_ani_coverage.txt > novel_ani.ini
cat Novel_ref.txt novel_ref_coverage.txt > novel_ref.ini


##Create the simulated metagenomes
conda activate /home/projects/ku-cbd/people/rapeis/.conda-anvio-7.1

cd ../../

../reads-for-assembly/gen-paired-end-reads 1_References/Novel_genomes/novel_ani.ini
../reads-for-assembly/gen-paired-end-reads 1_References/Novel_genomes/novel_ref.ini

gzip 2_Simulated_reads/Novel_genomes_*

# Run SingleM
cd 2_Simulated_reads

for i in Novel*-R1.fastq.gz;
do ../../singlem_versions/singlem_dev_03_08_22/bin/singlem pipe \
--singlem_metapackage ../0_Database/S3.metapackage_20220513.smpkg/ \
--forward $i \
--reverse ${i/-R1.fastq/-R2.fastq} \
--threads 40 \
--assignment-method naive_then_diamond \
--assignment-singlem-db ../0_Database/S3.metapackage_20220513.smpkg/gtdb_r207.reassigned.v5.sdb \
--archive-otu-table ../3_Outputs/5_Novel_genomes/${i/-R1.fastq.gz/_AUGUST_pipe.json} \
--otu-table ../3_Outputs/5_Novel_genomes/${i/-R1.fastq.gz/_AUGUST_pipe.tsv};
done


cd ../3_Outputs/5_Novel_genomes

for i in Novel*.json;
do ../../../singlem_versions/singlem_dev_03_08_22/bin/singlem condense \
--singlem_metapackage ../../0_Database/S3.metapackage_20220513.smpkg/ \
--input-archive-otu-table $i \
--output-otu-table ${i/.json/_condense.tsv} \
--apply-expectation-maximisation \
--trim-percent 10;
done
