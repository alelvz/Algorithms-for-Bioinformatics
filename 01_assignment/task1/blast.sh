#!/bin/bash
#SBATCH --job-name=blastn
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=30:00:00
#SBATCH --mem=32GB
#SBATCH --mail-type=fail

module load blast

# Create directories for databases and results
mkdir -p blast
mkdir -p results

# Database preparation
makeblastdb -in GCF_000005845.2_ASM584v2_genomic.fna -dbtype nucl -out blast/ecoli_db
makeblastdb -in GCF_000006765.1_ASM676v1_genomic.fna -dbtype nucl -out blast/paeruginosa_db
makeblastdb -in GCF_000009045.1_ASM904v1_genomic.fna -dbtype nucl -out blast/bsubtilis_db
makeblastdb -in GCF_000013425.1_ASM1342v1_genomic.fna -dbtype nucl -out blast/saureus_db
makeblastdb -in GCF_000195955.2_ASM19595v2_genomic.fna -dbtype nucl -out blast/mtuberculosis_db

# Reads and databases directories
read_dir="project1-data"
db_dir="blast"

# Run BLASTn for each genome and each read file
for genome in ecoli paeruginosa bsubtilis saureus mtuberculosis
do
    for read in miseq_10k_R1 miseq_10k_R2 no_errors_10k_R1 no_errors_10k_R2
    do
        echo "Aligning reads from $read to $genome"
        blastn -query ${read_dir}/${read}.fasta -db ${db_dir}/${genome}_db \
               -out results/${genome}_${read}_blastn.txt -outfmt 6 \
               -num_threads 16  -evalue 1e-6 -word_size 31

        echo "BLASTn alignment complete for $genome with $read"
    done
done

echo "All BLASTn alignments completed."

