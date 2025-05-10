#!/usr/bin/env python3
import argparse
from collections import defaultdict, Counter

def parse_fastq(fastq_file):
    sequences = []
    with open(fastq_file, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            sequences.append(seq)
    return sequences

def get_kmers(sequence, k):
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def build_debruijn_graph(kmers):
    edges = defaultdict(list)
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        if suffix not in edges[prefix]:
            edges[prefix].append(suffix)
            out_degree[prefix] += 1
            in_degree[suffix] += 1
    return edges, in_degree, out_degree

def get_largest_connected_component(edges):
    visited = set()
    largest_component = []
    def dfs_iterative(start):
        stack = [start]
        component = []
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                component.append(node)
                for neighbor in edges.get(node, []):
                    if neighbor not in visited:
                        stack.append(neighbor)
        return component
    all_nodes = set(edges.keys())
    for targets in edges.values():
        all_nodes.update(targets)
    for node in all_nodes:
        if node not in visited:
            component = dfs_iterative(node)
            if len(component) > len(largest_component):
                largest_component = component
    return largest_component

def hierholzer_algorithm(edges, start):
    edges_copy = {node: edges[node].copy() for node in edges}
    path = []
    stack = [start]
    while stack:
        current = stack[-1]
        if current in edges_copy and edges_copy[current]:
            stack.append(edges_copy[current].pop())
        else:
            path.append(stack.pop())
    path.reverse()
    return path if len(path) > 1 else []

def find_bubbles(edges, in_degree, out_degree):
    bubbles = []
    for node in edges:
        if len(edges[node]) > 1:
            visited = set()
            paths = []
            for next_node in edges[node]:
                path = [node, next_node]
                current = next_node
                seen = {node, next_node}
                while (
                    in_degree[current] == 1 and out_degree[current] == 1 and
                    current in edges and edges[current]
                ):
                    next_node = edges[current][0]
                    if next_node in seen:
                        break
                    path.append(next_node)
                    seen.add(next_node)
                    current = next_node
                if current not in visited:
                    paths.append(path)
                    visited.add(current)
            endpoints = defaultdict(list)
            for path in paths:
                endpoints[path[-1]].append(path)
            for end, ps in endpoints.items():
                if len(ps) > 1 and in_degree[end] > 1:
                    bubbles.extend(ps)
    return bubbles

def paths_to_contigs(paths, k):
    contigs = []
    for path in paths:
        if not path:
            continue
        contig = path[0] + ''.join(node[-1] for node in path[1:])
        contigs.append(contig)
    return contigs

def filter_kmers_by_coverage(kmers, min_coverage):
    counts = Counter(kmers)
    return [kmer for kmer, count in counts.items() if count >= min_coverage]

def write_fasta(contigs, output_file):
    with open(output_file, 'w') as f:
        for i, contig in enumerate(contigs):
            f.write(f">contig_{i+1}\n")
            for j in range(0, len(contig), 70):
                f.write(f"{contig[j:j+70]}\n")

def write_gfa(edges, output_file, bubbles=None):
    with open(output_file, 'w') as f:
        f.write("H\tVN:1.0\n")
        nodes = set(edges.keys())
        for targets in edges.values():
            nodes.update(targets)
        for node in nodes:
            f.write(f"S\t{node}\t*\n")
        for src, targets in edges.items():
            for tgt in targets:
                f.write(f"L\t{src}\t+\t{tgt}\t+\t0M\n")
        if bubbles:
            for i, path in enumerate(bubbles):
                f.write(f"# Bubble {i+1}: {'->'.join(path)}\n")

def main():
    parser = argparse.ArgumentParser(description='De Bruijn Graph Genome Assembler')
    parser.add_argument('-i', '--input', required=True, help='Input FASTQ file')
    parser.add_argument('-k', '--kmer', type=int, default=31, help='k-mer size')
    parser.add_argument('-o', '--output', default='contigs.fasta', help='Output FASTA file')
    parser.add_argument('-c', '--min-coverage', type=int, default=1, help='Min k-mer coverage')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    args = parser.parse_args()

    print(f"Reading sequences from {args.input}")
    sequences = parse_fastq(args.input)
    print(f"Found {len(sequences)} sequences")

    print(f"Extracting {args.kmer}-mers")
    all_kmers = [k for seq in sequences for k in get_kmers(seq, args.kmer)]
    print(f"Generated {len(all_kmers)} k-mers")

    kmers = filter_kmers_by_coverage(all_kmers, args.min_coverage) if args.min_coverage > 1 else all_kmers
    print(f"Kept {len(kmers)} k-mers after filtering")

    print("Building de Bruijn graph")
    edges, in_deg, out_deg = build_debruijn_graph(kmers)
    print(f"Graph has {len(edges)} nodes")

    print("Finding all connected components")
    def find_connected_components(edges):
        visited = set()
        components = []
        def dfs(start):
            stack = [start]
            comp = []
            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    comp.append(node)
                    for neighbor in edges.get(node, []):
                        if neighbor not in visited:
                            stack.append(neighbor)
            return comp
        all_nodes = set(edges.keys())
        for targets in edges.values():
            all_nodes.update(targets)
        for node in all_nodes:
            if node not in visited:
                components.append(dfs(node))
        return components

    paths = []
    components = find_connected_components(edges)
    print(f"Found {len(components)} components")

    for component in components:
        subgraph = {n: [x for x in edges.get(n, []) if x in component] for n in component if n in edges}
        sub_in = defaultdict(int)
        sub_out = defaultdict(int)
        for n, targets in subgraph.items():
            sub_out[n] = len(targets)
            for t in targets:
                sub_in[t] += 1
        start_node = next((n for n in component if sub_out[n] - sub_in[n] == 1), None)
        if not start_node and subgraph:
            start_node = next(iter(subgraph.keys()))
        if start_node:
            path = hierholzer_algorithm(subgraph, start_node)
            if path:
                paths.append(path)


    print("Detecting bubbles")
    bubble_paths = find_bubbles(edges, in_deg, out_deg)
    paths.extend(bubble_paths)

    print("Constructing contigs")
    contigs = paths_to_contigs(paths, args.kmer)
    contigs = list({c for c in contigs if len(c) >= args.kmer})
    print(f"Generated {len(contigs)} unique contigs")

    print(f"Writing contigs to {args.output}")
    write_fasta(contigs, args.output)

    gfa_out = args.output.rsplit('.', 1)[0] + '.gfa'
    print(f"Writing GFA graph to {gfa_out}")
    write_gfa(edges, gfa_out, bubbles=bubble_paths)

    print("Assembly complete!")

if __name__ == "__main__":
    main()
