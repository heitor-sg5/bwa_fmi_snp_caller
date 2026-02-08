from collections import defaultdict

from .aligner import extend_seed

def reverse_complement(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]

def generate_seeds(read, min_seed=12, max_seed=16):
    seeds = []
    seed_size = min(max_seed, len(read))
    mid = len(read) // 2
    start = max(0, mid - seed_size // 2)
    seeds.append((start, seed_size))
    fallback_seed_size = max(min_seed, seed_size // 2)
    step = max(1, fallback_seed_size // 2)
    for i in range(0, len(read) - fallback_seed_size + 1, step):
        seeds.append((i, fallback_seed_size))
    return seeds

def map_reads(reference, reads, fm, max_mismatches=3):
    alignments = defaultdict(list)
    for read_id, read in reads.items():
        best = None
        best_strand = "+"
        candidates = set()
        strands = {}
        for strand, seq in [("+", read), ("-", reverse_complement(read))]:
            seeds = generate_seeds(seq)
            for start, seed_len in seeds:
                seed = seq[start:start + seed_len]
                positions = fm.search_exact(seed)
                for pos in positions:
                    ref_pos = pos - start
                    if 0 <= ref_pos <= len(reference) - len(seq):
                        candidates.add(ref_pos)
                        strands[ref_pos] = strand
                if candidates:
                    break
        for p in candidates:
            strand = strands.get(p, "+")
            seq_to_align = read if strand == "+" else reverse_complement(read)
            aln = extend_seed(reference, seq_to_align, p, max_mismatches)
            if aln and (best is None or aln["mismatches"] < best["mismatches"]):
                best = aln
                best_strand = strand
        if best:
            best["strand"] = best_strand
            alignments[read_id].append(best)
    return alignments