from collections import defaultdict

from .aligner import extend_seed

def map_reads(reference, reads, fm, max_mismatches=3):
    alignments = defaultdict(list)
    for read_id, read in reads.items():
        candidates = set()
        seed_size = 12
        step = 4
        if len(read) >= seed_size:
            for i in range(0, len(read) - seed_size + 1, step):
                seed = read[i:i + seed_size]
                positions = fm.search_exact(seed)
                candidates.update(positions)
        
        if not candidates and len(read) >= 10:
            seed_size = 10
            step = 2
            for i in range(0, len(read) - seed_size + 1, step):
                seed = read[i:i + seed_size]
                positions = fm.search_exact(seed)
                candidates.update(positions)
        
        best = None
        for p in candidates:
            aln = extend_seed(reference, read, p, max_mismatches)
            if aln and (best is None or aln["mismatches"] < best["mismatches"]):
                best = aln
        
        if best:
            alignments[read_id].append(best)
    
    return alignments