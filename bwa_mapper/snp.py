from collections import defaultdict

def call_snps(alignments, reference, min_depth, min_alt_frac):
    pileup = defaultdict(lambda: defaultdict(int))
    for aln_list in alignments.values():
        for aln in aln_list:
            ref_pos = aln["ref_pos"]
            read_seq = aln["read_seq"]
            for i, base in enumerate(read_seq):
                pileup[ref_pos + i][base] += 1
    snps = []
    for pos in sorted(pileup.keys()):
        counts = pileup[pos]
        depth = sum(counts.values())
        if depth < min_depth:
            continue
        ref_base = reference[pos]
        for alt_base, cnt in counts.items():
            if alt_base != ref_base:
                alt_frac = cnt / depth
                if alt_frac >= min_alt_frac:
                    snps.append((pos, ref_base, alt_base, cnt, depth))
    return snps