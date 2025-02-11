import gzip
import sys

def read_fasta(file_path):
    """Reads a FASTA file and returns a dictionary with sequences per chromosome (converted to uppercase)."""
    sequences = {}
    current_chrom = None
    current_seq = []
    
    with gzip.open(file_path, "rt") if file_path.endswith(".gz") else open(file_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_chrom:
                    sequences[current_chrom] = "".join(current_seq).upper()  # Convert to uppercase
                current_chrom = line.strip().split()[0][1:]  # Remove ">" and take the chromosome name
                current_seq = []
            else:
                current_seq.append(line.strip())

    if current_chrom:
        sequences[current_chrom] = "".join(current_seq).upper()  # Convert to uppercase

    return sequences

def read_alu_sequence(file_path):
    """Reads the AluY sequence from a FASTA file and converts it to uppercase."""
    with open(file_path, "r") as f:
        lines = f.readlines()
    return "".join(line.strip() for line in lines[1:]).upper()  # Skip header and convert to uppercase

def exact_pattern_matching(genome, alu_seq):
    """Finds exact matches of AluY in a given genome sequence."""
    alu_len = len(alu_seq)
    matches = []

    for i in range(len(genome) - alu_len + 1):
        if genome[i:i+alu_len] == alu_seq:  # Character-by-character comparison (case insensitive)
            matches.append(i)

    return matches

def main(genome_file, alu_file):
    """Main function to run exact pattern matching."""
    print(f"Processing genome: {genome_file}")
    print(f"Using AluY sequence from: {alu_file}")

    alu_seq = read_alu_sequence(alu_file)  # Read and convert AluY to uppercase
    genome_sequences = read_fasta(genome_file)  # Read and convert genome to uppercase

    total_matches = 0

    for chrom, sequence in genome_sequences.items():
        print(f"Searching in {chrom}...")
        matches = exact_pattern_matching(sequence, alu_seq)
        total_matches += len(matches)
        
        for match in matches[:10]:  # Print first 10 matches per chromosome for sanity check
            print(f"Match at {chrom}: {match}-{match + len(alu_seq)}")

        print(f"Total matches in {chrom}: {len(matches)}\n")

    print(f"Total AluY matches in genome: {total_matches}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python exact_match.py <genome_fasta.gz> <alu_fasta>")
        sys.exit(1)

    genome_file = sys.argv[1]
    alu_file = sys.argv[2]

    main(genome_file, alu_file)
