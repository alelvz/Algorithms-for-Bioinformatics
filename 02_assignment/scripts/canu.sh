#!/bin/bash
#SBATCH --job-name=canu_ont
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00

module load canu

# Define input and output
INPUT_DIR=../data/synthetic_dataset/reads

# Error-free
#canu -p canu_noerror -d ../results/canu_noerror genomeSize=30k -nanopore-corrected ${INPUT_DIR}/no_error_ont_hq_50x.fastq maxThreads=16 maxMemory=64 useGrid=false -assemble stopOnLowCoverage=0 minInputCoverage=0

# With errors
#canu -p canu_error -d ../results/canu_error genomeSize=30k -nanopore ${INPUT_DIR}/ont_hq_50x.fastq maxThreads=16 maxMemory=64 useGrid=false

