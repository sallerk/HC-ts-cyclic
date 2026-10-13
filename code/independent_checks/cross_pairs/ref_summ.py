#!/usr/bin/env python3
# Summarise ref_xu.py outputs. Usage: python -I -B ref_summ.py file1 file2 ...
import sys, ast
from math import gcd
rows = []
for fn in sys.argv[1:]:
    with open(fn) as f:
        lines = f.read().splitlines()
    assert lines[-1] == "DONE", fn
    for L in lines[:-1]:
        rows.append(ast.literal_eval(L))
print("pairs (d,n):", len(rows), " d range", min(r["d"] for r in rows), max(r["d"] for r in rows))
bad = [r for r in rows if not (r["xu_ok"] and r["exact_ok"] and r["dim_ok"] and r.get("sigma_agrees", True)
                                and r.get("dim_ne_2k", True) and r.get("min_h11_minus_k", 0) >= 0
                                and r.get("z1fail", 0) == 0 and r.get("lift_missing", 0) == 0
                                and r.get("swaps_in_KR", True) and r.get("swapformula", True))]
print("rows with any failed check:", len(bad))
for r in bad[:20]:
    print("  BAD", r)
fails = sorted((r["d"], r["n"]) for r in rows if r["nfail"] > 0)
spor = [(d, n) for (d, n) in fails if n != d and 2 * n != d]
print("failing (d,n) other than n=d, n=d/2:", spor)
print("n=d/2 pairs with no failure:", [(r["d"], r["n"]) for r in rows if 2 * r["n"] == r["d"] and r["nfail"] == 0])
print("n=d with no failure:", [(r["d"], r["n"]) for r in rows if r["n"] == r["d"] and r["nfail"] == 0])
print("n|d-1 with failure:", [(r["d"], r["n"]) for r in rows if (r["d"] - 1) % r["n"] == 0 and r["nfail"] > 0])
mins = [r["min_h11_minus_k"] for r in rows if "min_h11_minus_k" in r]
print("Theorem C range: #pairs", len(mins), " min(h11-k) over all =", min(mins) if mins else None,
      " #pairs with min 0:", sum(1 for x in mins if x == 0))
kr = {r["d"]: (r["nfail"], r["kappa"]) for r in rows if r["n"] == r["d"]}
if kr:
    print("KR table d: |KR|(kappa)")
    print(" ".join("%d:%d(%d)" % (d, kr[d][0], kr[d][1]) for d in sorted(kr)))
    print("kappa=0 iff gcd(d,6)=1 :", all((kr[d][1] == 0) == (gcd(d, 6) == 1) for d in kr))
    def isprime(p):
        return p > 1 and all(p % q for q in range(2, int(p ** 0.5) + 1))
    print("d=2p, p>=11: kappa=48(p-1)?", [(d, kr[d][1], 48 * (d // 2 - 1)) for d in sorted(kr) if d % 2 == 0 and isprime(d // 2) and d // 2 >= 11 and kr[d][1] != 48 * (d // 2 - 1)])
    print("d=3p, p>=11: kappa=24(p-1)?", [(d, kr[d][1], 24 * (d // 3 - 1)) for d in sorted(kr) if d % 3 == 0 and isprime(d // 3) and d // 3 >= 11 and kr[d][1] != 24 * (d // 3 - 1)])
nd2 = [(r["d"], r["nfail"]) for r in rows if 2 * r["n"] == r["d"]]
print("n=d/2 failing ordered-pair counts (first 20):", nd2[:20])
print("sporadic failing counts:", [(r["d"], r["n"], r["nfail"]) for r in rows if r["nfail"] > 0 and r["n"] != r["d"] and 2 * r["n"] != r["d"]])
print("dims (6,2),(6,3),(6,6):", [(r["d"], r["n"], r["dimHdg"]) for r in rows if r["d"] == 6])
print("max sec:", max(r["sec"] for r in rows))
