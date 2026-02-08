from collections import defaultdict

from .aligner import extend_seed

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
        candidates = set()

        seeds = generate_seeds(read)

        for start, seed_len in seeds:
            seed = read[start:start + seed_len]
            positions = fm.search_exact(seed)
            for pos in positions:
                ref_pos = pos - start
                if 0 <= ref_pos <= len(reference) - len(read):
                    candidates.add(ref_pos)
            if candidates:
                break

        best = None
        for p in candidates:
            aln = extend_seed(reference, read, p, max_mismatches)
            if aln and (best is None or aln["mismatches"] < best["mismatches"]):
                best = aln

        if best:
            alignments[read_id].append(best)

    return alignments