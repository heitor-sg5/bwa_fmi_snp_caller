from collections import defaultdict

transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
alphabet = {"A", "C", "G", "T"}

def calculate_snp_stats(snp_list):
    total_snps = len(snp_list)
    ti_count = 0
    tv_count = 0
    alt_counts = defaultdict(int)

    for pos, ref, alt, count, depth in snp_list:
        alt_counts[alt] += count
        if (ref, alt) in transitions:
            ti_count += 1
        else:
            tv_count += 1
    ti_tv_ratio = (ti_count / tv_count) if tv_count > 0 else None
    most_common_alt = max(alt_counts.items(), key=lambda x: x[1])[0] if alt_counts else None
    stats = {
        "Total SNPs": total_snps,
        "Transitions": ti_count,
        "Transversions": tv_count,
        "Ti/Tv Ratio": round(ti_tv_ratio, 2) if ti_tv_ratio is not None else None,
        "Most Common Alt Alleles": most_common_alt
    }
    return stats

def calculate_mapping_stats(alignments, total_reads, reference_length):
    mapped_reads = sum(1 for aln_list in alignments.values() if aln_list)
    percent_mapped = (mapped_reads / total_reads) * 100 if total_reads > 0 else 0

    forward_reads = sum(1 for aln_list in alignments.values() for aln in aln_list if aln.get("strand") == "+")
    reverse_reads = sum(1 for aln_list in alignments.values() for aln in aln_list if aln.get("strand") == "-")
    forward_reverse_ratio = (forward_reads / reverse_reads) if reverse_reads > 0 else None

    stats = {
        "Reference Length": reference_length,
        "Total Reads": total_reads,
        "Mapped Reads": mapped_reads,
        "Percent Mapped": round(percent_mapped, 2),
        "Forward Strand Reads": forward_reads,
        "Reverse Strand Reads": reverse_reads,
        "Forward/Reverse Ratio": round(forward_reverse_ratio, 2) if forward_reverse_ratio is not None else None
    }
    return stats