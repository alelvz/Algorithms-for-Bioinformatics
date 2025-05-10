#!/usr/bin/env python3
"""
Overlap-Layout-Consensus (OLC) Assembler

This script implements an OLC assembler that:
1. Reads FASTQ files
2. Computes overlaps between reads (forward strands only)
3. Builds an overlap graph
4. Traverses the graph to find paths
5. Generates consensus sequences
6. Outputs contigs to a FASTA file
"""

import argparse
import os
import time
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Set, Optional


@dataclass
class Read:
    """Class for storing read information."""
    id: str
    sequence: str
    quality: str


@dataclass
class Overlap:
    """Class for storing overlap information between two reads."""
    read_a: str  # ID of first read
    read_b: str  # ID of second read
    overlap_len: int  # Length of overlap
    score: float  # Quality score for the overlap
    offset: int  # Offset of read_b relative to read_a
    mismatches: int  # Number of mismatches in the overlap


class OLCAssembler:
    def __init__(self, min_overlap: int = 20, max_mismatches: int = 2, min_contig_len: int = 200):
        """
        Initialize the OLC assembler.
        
        Args:
            min_overlap: Minimum overlap length between reads
            max_mismatches: Maximum number of mismatches allowed in overlaps
            min_contig_len: Minimum length for output contigs
        """
        self.min_overlap = min_overlap
        self.max_mismatches = max_mismatches
        self.min_contig_len = min_contig_len
        self.reads = {}  # Dict[str, Read]
        self.overlaps = []  # List[Overlap]
        self.graph = defaultdict(list)  # Dict[str, List[Tuple[str, Overlap]]]

    def parse_fastq(self, filename: str) -> None:
        """
        Parse a FASTQ file and store reads.
        
        Args:
            filename: Path to the FASTQ file
        """
        print(f"Parsing FASTQ file: {filename}")
        try:
            with open(filename, 'r') as file:
                lines = file.readlines()
                
            i = 0
            while i < len(lines):
                if i + 4 > len(lines):
                    break
                    
                header = lines[i].strip()
                if not header.startswith('@'):
                    i += 1
                    continue
                    
                read_id = header[1:].split()[0]
                sequence = lines[i + 1].strip()
                # Skip the '+' line
                quality = lines[i + 3].strip()
                
                self.reads[read_id] = Read(id=read_id, sequence=sequence, quality=quality)
                i += 4
                
            print(f"Parsed {len(self.reads)} reads")
        except Exception as e:
            print(f"Error parsing FASTQ file: {e}")
            exit(1)

    def compute_hamming_distance(self, seq1: str, seq2: str) -> int:
        """
        Compute the Hamming distance between two sequences of the same length.
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            
        Returns:
            Number of positions where the sequences differ
        """
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

    def find_overlap(self, seq1: str, seq2: str) -> Tuple[int, int, int]:
        """
        Find the best overlap between two sequences.
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            
        Returns:
            Tuple of (overlap length, number of mismatches, offset)
        """
        best_overlap = 0
        best_mismatches = float('inf')
        best_offset = 0
        
        # Try different overlap lengths, starting from the largest possible
        max_len = min(len(read_a_seq), len(read_b_seq))
        for overlap_len in range(max_len, self.min_overlap - 1, -1):

            # seq1 end overlaps with seq2 start
            suffix = seq1[-overlap_len:]
            prefix = seq2[:overlap_len]
            
            mismatches = self.compute_hamming_distance(suffix, prefix)
            
            if mismatches <= self.max_mismatches:
                best_overlap = overlap_len
                best_mismatches = mismatches
                best_offset = len(seq1) - overlap_len
                break  # Take the first good overlap we find (longest)
                
        return best_overlap, best_mismatches, best_offset

    def compute_overlaps(self) -> None:
        """
        Compute all-vs-all overlaps between reads.
        Performance optimized by using a suffix-prefix match approach.
        """
        print("Computing overlaps...")
        start_time = time.time()
        
        read_ids = list(self.reads.keys())
        read_count = len(read_ids)
        
        # Use a more efficient approach to report progress
        total_comparisons = read_count * (read_count - 1)
        progress_step = max(1, total_comparisons // 10)
        comparison_count = 0
        last_report_time = start_time
        
        # Index reads by their prefixes for faster overlap detection
        prefix_index = {}
        for read_id in read_ids:
            read = self.reads[read_id]
            # Store the minimum overlap length prefix for each read
            if len(read.sequence) >= self.min_overlap:
                prefix = read.sequence[:self.min_overlap]
                if prefix not in prefix_index:
                    prefix_index[prefix] = []
                prefix_index[prefix].append(read_id)
        
        # For each read, find potential overlaps
        for i, read_a_id in enumerate(read_ids):
            read_a = self.reads[read_a_id]
            read_a_seq = read_a.sequence
            
            # For each possible overlap length, check for matches
            for overlap_len in range(self.min_overlap, min(len(read_a_seq), 100) + 1):
                suffix = read_a_seq[-overlap_len:]
                
                # Check if any read has this suffix as a prefix
                potential_matches = prefix_index.get(suffix[:self.min_overlap], [])
                
                for read_b_id in potential_matches:
                    if read_a_id == read_b_id:
                        continue
                        
                    comparison_count += 1
                    
                    # Report progress at intervals
                    if comparison_count % progress_step == 0 or time.time() - last_report_time > 5:
                        progress = comparison_count / total_comparisons * 100
                        elapsed = time.time() - start_time
                        print(f"Progress: {progress:.1f}% - Elapsed: {elapsed:.1f}s - Found {len(self.overlaps)} overlaps")
                        last_report_time = time.time()
                    
                    read_b = self.reads[read_b_id]
                    read_b_seq = read_b.sequence
                    
                    # Verify the full overlap
                    if len(read_b_seq) >= overlap_len and read_b_seq[:overlap_len] == suffix:
                        # Perfect match
                        score = overlap_len
                        offset = len(read_a_seq) - overlap_len
                        
                        overlap = Overlap(
                            read_a=read_a_id,
                            read_b=read_b_id,
                            overlap_len=overlap_len,
                            score=score,
                            offset=offset,
                            mismatches=0
                        )
                        
                        self.overlaps.append(overlap)
                    else:
                        # Check with allowed mismatches
                        if len(read_b_seq) >= overlap_len:
                            mismatches = self.compute_hamming_distance(suffix, read_b_seq[:overlap_len])
                            
                            if mismatches <= self.max_mismatches:
                                # Calculate score (higher is better)
                                score = overlap_len / (mismatches + 1)
                                offset = len(read_a_seq) - overlap_len
                                
                                overlap = Overlap(
                                    read_a=read_a_id,
                                    read_b=read_b_id,
                                    overlap_len=overlap_len,
                                    score=score, 
                                    offset=offset,
                                    mismatches=mismatches
                                )
                                
                                self.overlaps.append(overlap)
        
        # Sort overlaps by score (higher score first)
        self.overlaps.sort(key=lambda x: (x.score, -x.mismatches), reverse=True)
        
        print(f"Found {len(self.overlaps)} overlaps in {time.time() - start_time:.1f} seconds")

    def build_overlap_graph(self) -> None:
        """
        Build a directed overlap graph from the computed overlaps.
        """
        print("Building overlap graph...")
        
        # Reset the graph
        self.graph = defaultdict(list)
        
        # Add edges for each overlap
        for overlap in self.overlaps:
            self.graph[overlap.read_a].append((overlap.read_b, overlap))
            
        # Count the number of nodes and edges
        num_nodes = len(set(node for node in self.graph.keys()) | 
                         set(edge for edges in self.graph.values() for edge, _ in edges))
        num_edges = sum(len(edges) for edges in self.graph.values())
        
        print(f"Built graph with {num_nodes} nodes and {num_edges} edges")

    def find_paths(self) -> List[List[Tuple[str, Optional[Overlap]]]]:
        """
        Find paths in the overlap graph using a greedy approach.
        
        Returns:
            List of paths, where each path is a list of (node, overlap) tuples
        """
        print("Finding paths in the overlap graph...")
        
        paths = []
        visited = set()
        
        # First, identify nodes with no incoming edges (potential starting points)
        all_nodes = set(self.graph.keys())
        for edges in self.graph.values():
            for next_node, _ in edges:
                all_nodes.add(next_node)
                
        incoming_edges = defaultdict(int)
        for node in self.graph:
            for next_node, _ in self.graph[node]:
                incoming_edges[next_node] += 1
                
        start_nodes = [node for node in all_nodes if incoming_edges[node] == 0]
        
        # If no clear starting nodes, use nodes with highest ratio of outgoing to incoming edges
        if not start_nodes:
            node_scores = {}
            for node in all_nodes:
                outgoing = len(self.graph[node]) if node in self.graph else 0
                incoming = incoming_edges[node]
                if incoming == 0:
                    node_scores[node] = float('inf')  # Prioritize nodes with no incoming edges
                else:
                    node_scores[node] = outgoing / incoming
                    
            sorted_nodes = sorted(node_scores.items(), key=lambda x: x[1], reverse=True)
            start_nodes = [node for node, _ in sorted_nodes[:min(len(sorted_nodes), 10)]]
        
        # Process each start node
        for start_node in start_nodes:
            if start_node in visited:
                continue
                
            path = [(start_node, None)]  # First node has no incoming overlap
            visited.add(start_node)
            
            current_node = start_node
            while True:
                # Find the best unvisited edge
                if current_node not in self.graph:
                    break
                    
                next_edges = [(next_node, overlap) for next_node, overlap in self.graph[current_node] 
                              if next_node not in visited]
                              
                if not next_edges:
                    break
                    
                # Sort by score
                next_edges.sort(key=lambda x: x[1].score, reverse=True)
                next_node, overlap = next_edges[0]
                
                path.append((next_node, overlap))
                visited.add(next_node)
                current_node = next_node
                
            if len(path) > 1:  # Only keep paths with at least 2 nodes
                paths.append(path)
        
        # Also find paths starting from any unvisited node
        remaining_nodes = all_nodes - visited
        for start_node in remaining_nodes:
            if start_node not in self.graph:
                continue
                
            path = [(start_node, None)]
            visited.add(start_node)
            
            current_node = start_node
            while True:
                # Find the best unvisited edge
                if current_node not in self.graph:
                    break
                    
                next_edges = [(next_node, overlap) for next_node, overlap in self.graph[current_node] 
                              if next_node not in visited]
                              
                if not next_edges:
                    break
                    
                # Sort by score
                next_edges.sort(key=lambda x: x[1].score, reverse=True)
                next_node, overlap = next_edges[0]
                
                path.append((next_node, overlap))
                visited.add(next_node)
                current_node = next_node
                
            if len(path) > 1:  # Only keep paths with at least 2 nodes
                paths.append(path)
                
        print(f"Found {len(paths)} paths")
        return paths

    def generate_consensus(self, path: List[Tuple[str, Optional[Overlap]]]) -> str:
        """
        Generate a consensus sequence from a path.
        
        Args:
            path: List of (node, overlap) tuples
            
        Returns:
            Consensus sequence
        """
        if not path:
            return ""
            
        # Start with the first read
        first_node, _ = path[0]
        first_read = self.reads[first_node]
        consensus = first_read.sequence
            
        # Add subsequent reads
        for i in range(1, len(path)):
            node, overlap = path[i]
            read = self.reads[node]
            
            # Extend the consensus by the non-overlapping part
            extension = read.sequence[overlap.overlap_len:]
            consensus += extension
            
        return consensus

    def assemble(self) -> List[str]:
        """
        Perform assembly by finding paths and generating consensus sequences.
        
        Returns:
            List of consensus sequences
        """
        paths = self.find_paths()
        
        print("Generating consensus sequences...")
        contigs = []
        
        for i, path in enumerate(paths):
            consensus = self.generate_consensus(path)
            
            # Only keep contigs that are longer than the minimum length
            if len(consensus) >= self.min_contig_len:
                contigs.append(consensus)
                
        # Sort contigs by length (longest first)
        contigs.sort(key=len, reverse=True)
                
        print(f"Generated {len(contigs)} contigs")
        return contigs

    def write_fasta(self, contigs: List[str], output_file: str) -> None:
        """
        Write contigs to a FASTA file.
        
        Args:
            contigs: List of contig sequences
            output_file: Path to the output FASTA file
        """
        print(f"Writing contigs to {output_file}")
        
        with open(output_file, 'w') as f:
            for i, contig in enumerate(contigs):
                f.write(f">contig_{i+1} length={len(contig)}\n")
                
                # Write sequence in lines of 60 characters
                for j in range(0, len(contig), 60):
                    f.write(contig[j:j+60] + "\n")
                    
        print(f"Wrote {len(contigs)} contigs to {output_file}")

    def run(self, input_file: str, output_file: str) -> None:
        """
        Run the full assembly pipeline.
        
        Args:
            input_file: Path to the input FASTQ file
            output_file: Path to the output FASTA file
        """
        start_time = time.time()
        
        print(f"Starting assembly with min_overlap={self.min_overlap}, max_mismatches={self.max_mismatches}")
        
        # Parse input FASTQ
        self.parse_fastq(input_file)
        
        # Compute overlaps
        self.compute_overlaps()
        
        # Build overlap graph
        self.build_overlap_graph()
        
        # Assemble contigs
        contigs = self.assemble()
        
        # Write output FASTA
        self.write_fasta(contigs, output_file)
        
        elapsed = time.time() - start_time
        print(f"Assembly completed in {elapsed:.2f} seconds")
        print(f"Generated {len(contigs)} contigs")
        
        # Print contig statistics
        if contigs:
            lengths = [len(contig) for contig in contigs]
            print(f"Longest contig: {max(lengths)} bp")
            print(f"Shortest contig: {min(lengths)} bp")
            print(f"Average contig length: {sum(lengths) / len(lengths):.1f} bp")
            print(f"Total assembly length: {sum(lengths)} bp")
        else:
            print("No contigs generated")


def main():
    """Main function to parse arguments and run the assembler."""
    parser = argparse.ArgumentParser(description='OLC Assembler')
    
    parser.add_argument('-i', '--input', required=True, help='Input FASTQ file')
    parser.add_argument('-o', '--output', default='contigs.fasta', help='Output FASTA file')
    parser.add_argument('-m', '--min-overlap', type=int, default=20, 
                        help='Minimum overlap length between reads')
    parser.add_argument('-M', '--max-mismatches', type=int, default=2,
                        help='Maximum number of mismatches allowed in overlaps')
    parser.add_argument('-c', '--min-contig-len', type=int, default=200,
                        help='Minimum length for output contigs')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.isfile(args.input):
        print(f"Error: Input file '{args.input}' does not exist")
        return
        
    # Run assembler
    assembler = OLCAssembler(
        min_overlap=args.min_overlap,
        max_mismatches=args.max_mismatches,
        min_contig_len=args.min_contig_len
    )
    
    assembler.run(args.input, args.output)


if __name__ == "__main__":
    main()