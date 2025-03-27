#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=5:00:00
#SBATCH --job-name=kraken_db
#SBATCH --mem=100Gb

# Download or prepare your reference genome files (FASTA format)
kraken2-build --download-taxonomy --db kraken2_custom_db
kraken2-build --add-to-library genomes/*.fasta --db kraken2_custom_db
kraken2-build --build --db kraken2_custom_db
kraken2-inspect --db kraken2_custom_db