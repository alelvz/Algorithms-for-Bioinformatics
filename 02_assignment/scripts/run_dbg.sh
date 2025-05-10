#!/bin/bash
#SBATCH --job-name=dbg_assembler
#SBATCH --output=logs/%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=30:00

K=40
outdir="../results/out_dbg/"
module load Bandage

# Process toy_dataset
for fq in ../data/toy_dataset/*.fastq; do
    base=$(basename "$fq" .fastq)
    echo "[✓] Processing $base"
    python dbg.py -i "$fq" -k $K -o "$outdir/${base}_dbg_k${K}.fasta"
    Bandage image $outdir/${base}_dbg_k${K}.gfa $outdir/${base}_dbg_k${K}.png
done

# Process synthetic_dataset/reads
for fq in ../data/synthetic_dataset/reads/*.fastq; do
    base=$(basename "$fq" .fastq)
    echo "[✓] Processing $base"
    python dbg.py -i "$fq" -k $K -o "$outdir/${base}_dbg_k${K}.fasta" -c 2
    Bandage image $outdir/${base}_dbg_k${K}.gfa $outdir/${base}_dbg_k${K}.png
done
