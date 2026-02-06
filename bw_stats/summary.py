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