import os
import gzip
import psutil
from collections import defaultdict, Counter

def load_genome(file_path):
    """ Load a genome sequence from a gzip-compressed FASTA file. """
    with gzip.open(file_path, 'rt') as file:
        # Skip the header line
        file.readline()
        # Read and concatenate the rest of the file, removing newline characters
        return file.read().replace('\n', '')

def build_kmer_index(genomes_dir, k):
    kmer_index = defaultdict(list)
    for filename in os.listdir(genomes_dir):
        if filename.endswith(".fna.gz"):
            file_path = os.path.join(genomes_dir, filename)
            sequence = load_genome(file_path)
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                kmer_index[kmer].append(filename)
    return kmer_index


def classify_reads(reads_file, kmer_index, k):
    """ Classify reads by searching each read's k-mers in the k-mer index. 
        Count the number of matches for each organism per read.
    """
    # This dictionary will hold the count of k-mers matching each genome.
    organism_kmer_count = defaultdict(Counter)
    
    with gzip.open(reads_file, 'rt') as file:
        for line in file:
            if line.startswith('@'):  # Skip header
                read_id = next(file).strip()  # Read the sequence line
                read_kmer_counts = Counter()  # Track k-mers found in this read
                for i in range(len(read_id) - k + 1):
                    kmer = read_id[i:i+k]
                    if kmer in kmer_index:
                        # Update k-mer counts for organisms
                        for organism in kmer_index[kmer]:
                            read_kmer_counts[organism] += 1
                # Update global counts for each organism
                for organism, count in read_kmer_counts.items():
                    organism_kmer_count[organism][read_id] += count
                next(file)  # Skip the '+' line
                next(file)  # Skip the quality score line
    return organism_kmer_count

k=31
genomes_dir = 'genomes'
process = psutil.Process()
start_memory = process.memory_info().rss 

kmer_index = build_kmer_index(genomes_dir, k)
indexed_memory = process.memory_info().rss 

data_dir = 'project1-data'

num_kmers = len(kmer_index)
print(f"There are {num_kmers} unique k-mers in the index.")

output_file = 'kmer_classification.out'

# Process each reads file in the directory
with open(output_file, 'w') as out_file:
    for reads_filename in os.listdir(data_dir):
        if reads_filename.endswith('.fastq.gz'):
            reads_file = os.path.join(data_dir, reads_filename)
            classified_results = classify_reads(reads_file, kmer_index, k)
            out_file.write(f'Results for {reads_filename}:\n')
            for organism, reads in classified_results.items():
                total_kmers = sum(reads.values())
                out_file.write(f"{organism} has {len(reads)} reads with k-mer matches.\n")
                out_file.write(f"Total k-mer matches for {organism}: {total_kmers}\n")
            out_file.write('\n')

classified_memory = process.memory_info().rss  # Memory after classification
memory_used = (classified_memory - start_memory) / (1024 ** 2)  # Convert to MB

print(f"Classification completed. Results are saved in {output_file}.")
print(f"Memory used: {memory_used:.2f} MB")