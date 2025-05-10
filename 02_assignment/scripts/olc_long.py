#!/usr/bin/env python3
"""
Optimized Overlap-Layout-Consensus (OLC) Assembly for Oxford Nanopore reads.
Includes prefix-based filtering to reduce redundant overlap comparisons.
"""

import argparse
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import networkx as nx
from Bio import SeqIO

class OLCAssembler:
    def __init__(self, min_overlap=500, min_identity=0.80, threads=4):
        self.min_overlap = min_overlap
        self.min_identity = min_identity
        self.threads = threads
        self.reads = {}
        self.overlap_graph = nx.DiGraph()

    def load_fastq(self, fastq_file):
        print(f"Loading reads from {fastq_file}...")
        count = 0
        for record in SeqIO.parse(fastq_file, "fastq"):
            self.reads[record.id] = str(record.seq)
            count += 1
        print(f"Loaded {count} reads")
        return count

    def generate_candidate_pairs(self, k=6):
        """Reduce comparison space using shared k-mer prefixes"""
        index = defaultdict(list)
        for read_id, seq in self.reads.items():
            if len(seq) >= k:
                index[seq[:k]].append(read_id)

        candidate_pairs = set()
        for group in index.values():
            for i in range(len(group)):
                for j in range(len(group)):
                    if i != j:
                        candidate_pairs.add((group[i], group[j]))
        return list(candidate_pairs)

    def compute_overlap(self, read_pair):
        read1_id, read2_id = read_pair
        if read1_id == read2_id:
            return None

        read1 = self.reads[read1_id]
        read2 = self.reads[read2_id]
        best_overlap = None
        best_score = 0

        for olen in range(min(len(read1), len(read2)), self.min_overlap - 1, -1):
            suffix = read1[-olen:]
            prefix = read2[:olen]
            matches = sum(s == p for s, p in zip(suffix, prefix))
            identity = matches / olen
            if identity >= self.min_identity and identity > best_score:
                best_score = identity
                best_overlap = (read1_id, read2_id, olen, identity)
                if identity > 0.95:
                    break
        return best_overlap

    def build_overlap_graph(self):
        print("Building overlap graph...")

        pairs = self.generate_candidate_pairs(k=6)

        print(f"[✓] Filtered candidate pairs: {len(pairs)}")

        overlaps = []
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            for result in executor.map(self.compute_overlap, pairs, chunksize=500):
                if result:
                    overlaps.append(result)

        for r1, r2, olen, score in overlaps:
            self.overlap_graph.add_edge(r1, r2, overlap=olen, score=score)
        print(f"Graph has {self.overlap_graph.number_of_edges()} edges")

    def find_paths(self):
        print("Finding non-branching paths...")
        paths = []
        visited = set()

        for node in self.overlap_graph.nodes():
            if node in visited or self.overlap_graph.in_degree(node) == 1:
                continue
            path = [node]
            visited.add(node)
            current = node
            while self.overlap_graph.out_degree(current) == 1:
                succ = list(self.overlap_graph.successors(current))[0]
                if succ in visited or self.overlap_graph.in_degree(succ) != 1:
                    break
                path.append(succ)
                visited.add(succ)
                current = succ
            if len(path) > 1:
                paths.append(path)

        for comp in nx.weakly_connected_components(self.overlap_graph):
            sub = self.overlap_graph.subgraph(comp)
            if all(sub.in_degree(n) == 1 and sub.out_degree(n) == 1 for n in sub.nodes):
                start = list(comp)[0]
                if start not in visited:
                    path = [start]
                    visited.add(start)
                    curr = list(sub.successors(start))[0]
                    while curr != start:
                        path.append(curr)
                        visited.add(curr)
                        curr = list(sub.successors(curr))[0]
                    paths.append(path)

        print(f"Found {len(paths)} paths")
        return paths

    def generate_layout(self, paths):
        print("Generating layouts...")
        layouts = []
        for path in paths:
            layout = []
            for i in range(len(path) - 1):
                r1 = path[i]
                r2 = path[i + 1]
                olen = self.overlap_graph[r1][r2]['overlap']
                if i == 0:
                    layout.append((r1, 0, len(self.reads[r1])))
                layout.append((r2, olen, len(self.reads[r2]) - olen))
            layouts.append(layout)
        return layouts

    def compute_consensus(self, layouts):
        print("Computing consensus sequences...")
        contigs = []
        for i, layout in enumerate(layouts):
            consensus = ""
            for j, (read_id, offset, length) in enumerate(layout):
                seq = self.reads[read_id]
                if j == 0:
                    consensus += seq
                else:
                    consensus += seq[offset:]
            contigs.append((f"contig_{i+1}", consensus))
        print(f"Generated {len(contigs)} contigs")
        return contigs

    def write_fasta(self, contigs, output_file):
        print(f"Writing to {output_file}...")
        with open(output_file, 'w') as f:
            for cid, seq in contigs:
                f.write(f">{cid} length={len(seq)}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")
        print(f"Wrote {len(contigs)} contigs")

    def assemble(self, fastq_file, output_file):
        self.load_fastq(fastq_file)
        self.build_overlap_graph()
        paths = self.find_paths()
        layouts = self.generate_layout(paths)
        contigs = self.compute_consensus(layouts)
        self.write_fasta(contigs, output_file)
        return len(contigs)


def main():
    parser = argparse.ArgumentParser(description='OLC Assembly for ONT reads')
    parser.add_argument('input', help='Input FASTQ file')
    parser.add_argument('output', help='Output FASTA file')
    parser.add_argument('-m', '--min-overlap', type=int, default=500,
                        help='Minimum overlap length (default: 500)')
    parser.add_argument('-i', '--min-identity', type=float, default=0.8,
                        help='Minimum identity in overlaps (default: 0.8)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='Number of CPU threads to use (default: 4)')
    args = parser.parse_args()

    assembler = OLCAssembler(
        min_overlap=args.min_overlap,
        min_identity=args.min_identity,
        threads=args.threads
    )
    assembler.assemble(args.input, args.output)

if __name__ == "__main__":
    main()
