#!/bin/bash
#SBATCH --job-name=kraken2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --mem=300GB
#SBATCH --mail-type=fail

db_dir='/ibex/scratch/projects/c2014/rund/metagenomics_prelimdata/pilot/output'

kraken2 --db ${db_dir} --report kraken_assemblies/${1}_kraken_report.txt --output kraken_assemblies/${1}_kraken_out.txt ../02_megahit/${1}/final.contigs.fa 
bracken -d ${db_dir} -i kraken_assemblies/${1}_kraken_report.txt -o bracken/${1}_bracken_out.txt -w bracken/${1}_bracken_report.txt -r 300 -l P

#kraken2 --db ${db_dir} --output test_kraken.txt --report test_report.txt ../02_megahit/${1}/final.contigs.fa 
#bracken -d ${db_dir} -i test_report.txt -o test_bracken.txt -w test_bracken.report -r 300 -l P 

