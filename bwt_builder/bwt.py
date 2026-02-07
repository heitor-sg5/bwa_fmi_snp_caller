def build_suffix_array(text):
    n = len(text)
    if n == 0:
        return []
    rank = [ord(c) for c in text]
    tmp = [0] * n
    sa = list(range(n))

    k = 1
    while k < n:
        def key(i):
            r2 = rank[i + k] if i + k < n else -1
            return (rank[i], r2)
        sa.sort(key=key)

        tmp[sa[0]] = 0
        for i in range(1, n):
            prev, curr = sa[i - 1], sa[i]
            prev_key = (rank[prev], rank[prev + k] if prev + k < n else -1)
            curr_key = (rank[curr], rank[curr + k] if curr + k < n else -1)
            tmp[curr] = tmp[prev] + (1 if prev_key != curr_key else 0)
        rank, tmp = tmp, rank
        k *= 2
        if rank[sa[-1]] == n - 1:
            break
    return sa

def build_bwt(text, sa, sentinel="$"):
    if sentinel not in text:
        text = text + sentinel
    bwt_chars = []
    for idx in sa:
        if idx == 0:
            bwt_chars.append(sentinel)
        else:
            bwt_chars.append(text[idx - 1])
    return "".join(bwt_chars)