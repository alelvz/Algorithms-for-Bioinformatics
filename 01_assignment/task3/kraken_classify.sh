#!/bin/bash
#SBATCH --job-name=kraken2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --mem=16GB
#SBATCH --mail-type=fail


module load kraken2

db_dir='/ibex/scratch/projects/c2014/alelopezv/2502_Algorithms/01_assignment/task3/kraken2_custom_db'
reads_dir='../project1-data'

mkdir -p kraken_results
mkdir -p bracken_results

# Loop through all R1 files and find corresponding R2
for R1 in ${reads_dir}/*_R1.fastq.gz; do
    # Infer R2 file name
    R2=${R1/_R1/_R2}
    
    # Extract sample name
    sample=$(basename ${R1} 10k_R1.fastq.gz)
    sample=${sample#simulated_reads_}

    # Run Kraken2 for paired-end reads
    kraken2 --db ${db_dir} \
            --paired \
            --report ${sample}_kraken2_report.txt \
            --output ${sample}_kraken2_output.txt \
            ${R1} ${R2}
done

# Loop through all paired-end R1 reads
for R1 in ${reads_dir}/*_R1.fastq.gz; do
    R2=${R1/_R1/_R2}
    
    # Clean sample name
    sample=$(basename ${R1} _10k_R1.fastq.gz)
    sample=${sample#simulated_reads_}

    echo "Processing ${sample}..."

    # Run Kraken2
    kraken2 --db ${db_dir} \
            --paired \
            --report kraken_results/${sample}_kraken_report.txt \
            --output kraken_results/${sample}_kraken_out.txt \
            ${R1} ${R2}

done