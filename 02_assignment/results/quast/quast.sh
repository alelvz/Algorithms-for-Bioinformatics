#!/bin/bash
#SBATCH --job-name=quast_eval
#SBATCH --output=%j_quast.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

module load quast

# Paths to reference genomes
REF_TOY_R="../../data/toy_dataset/reference_r.fasta"
REF_TOY_B="../../data/toy_dataset/reference_b.fasta"
REF_MERS="../../data/synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna"

mkdir -p quast_results

# Function to run QUAST with a given reference
run_quast() {
    local asm=$1
    local ref=$2
    local label=$(basename "$asm" .fasta)
    echo "[✓] Evaluating $label against $(basename $ref)"
    quast.py -r "$ref" -o "${label}" "$asm" --threads 8 --labels "$label" -m 300
}

# Analyze DBG assemblies
for f in ../out_dbg/*.fasta; do
    fname=$(basename "$f")
    if [[ "$fname" == *r_dbg* ]]; then
        run_quast "$f" "$REF_TOY_R"
    elif [[ "$fname" == *b_dbg* ]]; then
        run_quast "$f" "$REF_TOY_B"
    else
        run_quast "$f" "$REF_MERS"
    fi
done

# Analyze OLC assemblies
for f in ../out_olc/*.fasta; do
    fname=$(basename "$f")
    if [[ "$fname" == *r_olc* ]]; then
        run_quast "$f" "$REF_TOY_R"
    elif [[ "$fname" == *b_olc* ]]; then
        run_quast "$f" "$REF_TOY_B"
    else
        run_quast "$f" "$REF_MERS"
    fi
done

