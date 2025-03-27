import sys
import gzip
import ahocorasick
from Bio import SeqIO
from collections import defaultdict
import time
import psutil

def load_genome(file_path):
    """Load genome sequence from a compressed FASTA file."""
    with gzip.open(file_path, "rt") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        return str(records[0].seq)

def build_trie(genome_files, k):
    """Build an Aho-Corasick trie from reference genomes."""
    A = ahocorasick.Automaton()
    for genome_id, file_path in genome_files.items():
        sequence = load_genome(file_path)
        for i in range(len(sequence) - k + 1):
            k_mer = sequence[i:i+k]
            if k_mer in A:
                existing_value = A.get(k_mer)
                updated_genome_ids = existing_value[1]
                updated_genome_ids.add(genome_id)
                A.remove_word(k_mer)
                A.add_word(k_mer, (k_mer, updated_genome_ids))
            else:
                A.add_word(k_mer, (k_mer, {genome_id}))
    A.make_automaton()
    return A

def classify_reads(reads_file, trie, k):
    """Classify reads using the trie and handle multiple matches."""
    read_matches = defaultdict(int)
    with gzip.open(reads_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            read_seq = str(record.seq)
            matches = set()
            for end_idx, (k_mer, genome_ids) in trie.iter(read_seq):
                matches.update(genome_ids)
            for match in matches:
                read_matches[match] += 1
    return read_matches

def main():
    # Start timing and memory measurement
    start_time = time.time()
    process = psutil.Process()
    
    if len(sys.argv) < 2:
        print("Usage: python script.py <file_suffix>")
        sys.exit(1)

    file_suffix = sys.argv[1]

    # Parameters
    k = 31  # k-mer size

    # Load reference genomes
    genome_files = {
        'E_coli': 'GCF_000005845.2_ASM584v2_genomic.fna.gz',
        'B_subtilis': 'GCF_000009045.1_ASM904v1_genomic.fna.gz',
        'P_aeruginosa': 'GCF_000006765.1_ASM676v1_genomic.fna.gz',
        'S_aureus': 'GCF_000013425.1_ASM1342v1_genomic.fna.gz',
        'M_tuberculosis': 'GCF_000195955.2_ASM19595v2_genomic.fna.gz'
    }

    # Build the trie from reference genomes
    trie = build_trie(genome_files, k)

    # Classify reads from one of the read files (can loop or modify for others)
    read_file = f"project1-data/simulated_reads_{file_suffix}.fastq.gz"
    classified_reads = classify_reads(read_file, trie, k)

    # Stop timing and memory measurement
    end_time = time.time()
    memory_usage = process.memory_info().rss  # resident set size; real memory usage

    # Print results (or handle otherwise)
    for organism, count in classified_reads.items():
        print(f"Reads matching {organism}: {count}")

    print(f"Execution time: {end_time - start_time:.2f} seconds")
    print(f"Memory usage: {memory_usage / (1024 * 1024):.2f} MB")  # Convert bytes to MB

if __name__ == "__main__":
    main()
