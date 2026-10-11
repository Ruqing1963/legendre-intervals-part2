"""
stable_stats.py -- statistics of the Frobenius-stable translation classes of hexads (data/hexads_k{k}.json):
for each S_6-class (type of the canonical sextic B), the distribution of the number of stable classes and the multiset
of factorisation types of the normalised sextics of all stable classes.  Also prints, for k = 3, the table of all pairs
(l3, l5) with their class and W' for the appendix of Part III.
Usage: python stable_stats.py 3 5
"""
import json
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def run(k):
    with open(os.path.join(HERE, "..", "data", f"hexads_k{k}.json")) as f:
        data = json.load(f)
    if isinstance(data, dict):
        items = list(data.items())
    else:
        items = [(None, d) for d in data]
    print(f"\nk={k}: {len(items)} records; keys of first record:", sorted(items[0][1].keys()) if isinstance(items[0][1], dict) else type(items[0][1]))
    import legendre_d6_qf as QF
    F = QF.Fq(1 << k)
    M, INV = F.M, F.INV
    nst = defaultdict(Counter)
    types = defaultdict(Counter)
    canon = Counter()
    for key, rec in items:
        if not isinstance(rec, dict) or "stable" not in rec:
            continue
        st = rec["stable"]
        l3, l5 = rec["l3"], rec["l5"]
        i5 = int(INV[l5])
        B = [int(M[i5, i5]), i5, int(M[int(M[l3, l3]), int(M[i5, i5])]), 0]      # b0, b1, b2, b3 of B/l5^2
        can = [tuple(s["sextic_type"]) for s in st if list(s["sextic"]) == B]
        assert len(can) == 1, (l3, l5, B, [s["sextic"] for s in st])
        label = can[0]
        nst[label][len(st)] += 1
        types[label][tuple(sorted(tuple(s["sextic_type"]) for s in st))] += 1
        canon[label] += 1
    for label in sorted(nst, key=str):
        print(f"  class {label}: #pairs {canon[label]}; #stable classes -> #pairs {dict(nst[label])}")
        for t, n in sorted(types[label].items(), key=str)[:6]:
            print(f"      stable sextic types {t}  x{n}")


if __name__ == "__main__":
    for k in [int(a) for a in sys.argv[1:]] or [3, 5]:
        run(k)
