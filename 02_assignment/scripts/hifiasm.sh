#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=12
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --account=cs249

module load hifiasm
dir='/ibex/reference/course/cs249/lizard/input'

# Run Hifiasm with PacBio HiFi + Hi-C
hifiasm -o lizard_asm_hic -t 48 \
  --h1 ${dir}/hic/lizard_hic_R1.fastq.gz \
  --h2 ${dir}/hic/lizard_hic_R2.fastq.gz \
  ${dir}/pacbio/lizard_liver.fastq.gz

