def hamming_distance(a, b):
    return sum(x != y for x, y in zip(a, b))

def extend_seed(reference, read, ref_pos, max_mismatches):
    if ref_pos + len(read) > len(reference):
        return None
    ref_segment = reference[ref_pos:ref_pos + len(read)]
    mismatches = hamming_distance(ref_segment, read)
    if mismatches <= max_mismatches:
        return {
            "ref_pos": ref_pos,
            "mismatches": mismatches,
            "read_seq": read
        }
    return None