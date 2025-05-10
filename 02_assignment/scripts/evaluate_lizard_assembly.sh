#!/bin/bash
#SBATCH --job-name=lizard_eval
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00

# Load modules
module load quast
module load Busco
module load merqury
module load flagger
module load gfatools

# Set input/output paths
ASSEMBLY_GFA="../results/hifiasm/lizard_asm_hic.hic.p_ctg.gfa"
ASSEMBLY_FASTA="../results/hifiasm/lizard_p_ctg.fasta"
REF_KMER_DB="../results/hifiasm_evaluation/merqury/read_hifi.meryl" 
READS='/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver.fastq.gz'


# Convert GFA to FASTA
gfatools gfa2fa $ASSEMBLY_GFA > $ASSEMBLY_FASTA

# QUAST
quast.py $ASSEMBLY_FASTA -o ../results/hifiasm_evaluation/quast -t 16

# BUSCO (for completeness evaluation)
busco -i ../results/hifiasm/lizard_p_ctg.fasta -o ../results/hifiasm_evaluation/busco -l sauropsida_odb10 -m genome --cpu 16

# Merqury (for QV and k-mer completeness)
meryl count k=21 output $REF_KMER_DB $READS
merqury.sh $REF_KMER_DB $ASSEMBLY_FASTA ../results/hifiasm_evaluation/merqury/

# Flagger (for misassembly detection)
flagger -r $ASSEMBLY_FASTA -o flagger_out --threads 16

