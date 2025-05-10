#!/bin/bash
#SBATCH --job-name=flye_ont
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00

module load flye

# Input and output
READS=../data/synthetic_dataset/reads
OUTDIR=../results/

# Run Flye for ONT data (no error)
#flye --nano-raw $READS/no_error_ont_hq_50x.fastq --genome-size 30k --out-dir $OUTDIR/flye_noerror --threads 16

# Run Flye for ONT data (error)
flye --nano-raw $READS/ont_hq_50x.fastq --genome-size 30k --out-dir $OUTDIR/flye_error --threads 16