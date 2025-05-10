#!/bin/bash
#SBATCH --job-name=quast_eval
#SBATCH --output=%j_quast.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

module load quast

REF_MERS="../../data/synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna"

# QUAST evaluations for MERS assemblies
quast.py ../spades_noerror/spades_noerror.fasta \
    -r $REF_MERS \
    -o spades_noerror \
    -t 8

quast.py ../spades_error/spades_error.fasta \
    -r $REF_MERS \
    -o spades_error \
    -t 8

quast.py ../flye_noerror/flye_noerror.fasta \
    -r $REF_MERS \
    -o flye_noerror \
    -t 8

quast.py ../flye_error/flye_error.fasta \
    -r $REF_MERS \
    -o flye_error \
    -t 8