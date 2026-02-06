from collections import defaultdict

def call_snps(alignments, reference, min_depth=3, min_alt_frac=0.2):
    pileup = defaultdict(lambda: defaultdict(int))
    for aln_list in alignments.values():
        for aln in aln_list:
            ref_pos = aln["ref_pos"]
            read_seq = aln["read_seq"]
            ref_seq = reference[ref_pos:ref_pos + len(read_seq)]
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

def format_snps(reference, snps):
    snp_dict = defaultdict(lambda: defaultdict(int))
    for pos, ref, alt, count in snps:
        snp_dict[pos][ref] = 0
        snp_dict[pos][alt] = count
    snp_list = []
    for pos in sorted(snp_dict.keys()):
        alleles = snp_dict[pos]
        ref_base = reference[pos]
        for alt_base, count in alleles.items():
            if alt_base != ref_base:
                snp_list.append((pos, ref_base, alt_base, count))
    return snp_list