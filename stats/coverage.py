def compute_coverage(alignments, ref_length):
    coverage = [0] * ref_length
    for aln_list in alignments.values():
        for aln in aln_list:
            start = aln["ref_pos"]
            end = start + len(aln["read_seq"])
            for i in range(start, min(end, ref_length)):
                coverage[i] += 1
    return coverage