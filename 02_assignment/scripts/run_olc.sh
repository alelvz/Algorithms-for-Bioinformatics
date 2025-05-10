#!/bin/bash
#SBATCH --job-name=olc_assembler
#SBATCH --output=logs/%j.out
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=30:00


outdir="out_olc"
mkdir -p $outdir

# Process toy_dataset
for fq in data/toy_dataset/*.fastq; do
    base=$(basename "$fq" .fastq)
    echo "[✓] Processing $base"
    python olc.py -i "$fq" -o "$outdir/${base}_olc.fasta" -m 30 -M 2 -c 150
done

# Process synthetic_dataset/reads
for fq in data/synthetic_dataset/reads/*hiseq*; do
    base=$(basename "$fq" .fastq)
    echo "[✓] Processing $base"
    python olc.py -i "$fq" -o "$outdir/${base}_olc.fasta" -m 40 -M 3 -c 200
done

# Process synthetic_dataset/reads
for fq in data/synthetic_dataset/reads/*ont*; do
    base=$(basename "$fq" .fastq)
    echo "[✓] Processing $base"
    python olc_long.py "$fq" "$outdir/${base}_olc.fasta" -m 200 -i 0.8
done
