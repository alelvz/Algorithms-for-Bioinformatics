import os
import gzip
import psutil
from collections import defaultdict, Counter

def load_genome(file_path):
    """ Load a genome sequence from a gzip-compressed FASTA file. """
    with gzip.open(file_path, 'rt') as file:
        file.readline()  # Skip the header line
        return file.read().replace('\n', '')
    
def get_minimizers(sequence, k, w):
    minimizers = set()
    for i in range(len(sequence) - k - w + 2):
        window = sequence[i:i + k + w - 1]
        kmers = [window[j:j + k] for j in range(w)]
        minimizer = min(kmers)
        minimizers.add(minimizer)
    return minimizers

def build_minimizer_index(genomes_dir, k, w):
    index = defaultdict(list)
    for filename in os.listdir(genomes_dir):
        if filename.endswith(".fna.gz"):
            file_path = os.path.join(genomes_dir, filename)
            sequence = load_genome(file_path)
            for minimizer in get_minimizers(sequence, k, w):
                index[minimizer].append(filename)
    return index

def classify_reads_minimizer(reads_file, minimizer_index, k, w):
    organism_counts = defaultdict(Counter)
    with gzip.open(reads_file, 'rt') as file:
        for line in file:
            if line.startswith('@'):
                read = next(file).strip()
                read_minimizers = get_minimizers(read, k, w)
                for minimizer in read_minimizers:
                    if minimizer in minimizer_index:
                        for org in minimizer_index[minimizer]:
                            organism_counts[org][read] += 1
                next(file)  # Skip '+'
                next(file)  # Skip quality
    return organism_counts

# Define k, w, and directories
k = 31
w = 10  # Minimizer window size
genomes_dir = 'genomes'
data_dir = 'project1-data'
output_file = 'kmer_minimizer.out'

# Memory tracking
process = psutil.Process()
start_mem = process.memory_info().rss

# Build index
min_index = build_minimizer_index(genomes_dir, k, w)
mid_mem = process.memory_info().rss

# Classify
with open(output_file, 'w') as out_file:
    for reads_filename in os.listdir(data_dir):
        if reads_filename.endswith('.fastq.gz'):
            reads_path = os.path.join(data_dir, reads_filename)
            result = classify_reads_minimizer(reads_path, min_index, k, w)
            out_file.write(f"Results for {reads_filename}:\n")
            for organism, reads in result.items():
                total = sum(reads.values())
                out_file.write(f"{organism} has {len(reads)} reads with k-mer matches.\n")
                out_file.write(f"Total k-mer matches for {organism}: {total}\n")
            out_file.write('\n')

end_mem = process.memory_info().rss
mem_usage = (end_mem - start_mem) / (1024**2)
print(f"Results saved to {output_file}")
print(f"Memory used: {mem_usage:.2f} MB")

