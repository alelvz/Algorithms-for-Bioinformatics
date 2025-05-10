#!/bin/bash
#SBATCH --job-name=spades
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=spades_%j.out
#SBATCH --error=spades_%j.err

module load spades
INPUT_DIR=../data/synthetic_dataset/reads

# Error-free
spades.py -s ${INPUT_DIR}/no_error_reads_hiseq_5k.fastq -o ../results/spades_noerror -t 16 -m 32 --only-assembler

# With errors
spades.py -s ${INPUT_DIR}/reads_hiseq_5k.fastq -o  ../results/spades_error -t 16 -m 32

