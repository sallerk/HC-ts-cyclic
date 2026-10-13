#!/usr/bin/env python3
# Referee's independent check for THEOREM YY.
# Uses ONLY the definitions of theorem_u.md: c_n (Sec. 2 Notation, as a direct count), Prop. 2.3 (h20, h02, dim),
# 0.3 (Galois action of (Z/L)^*), 0.4 (pure / non-pure / zero, pure orbit), Theorem U(4) for (XU).
# None of the manuscript's formulas (N*, sigma, Lemma 2.5) is used to decide (XU); sigma is computed only to compare.
# Usage: python -I -B ref_xu.py full DMIN DMAX OUTFILE [mlist]
#        python -I -B ref_xu.py fam DMAX OUTFILE
import sys, time
import numpy as np
from math import gcd

def lcm(a, b):
    return a * b // gcd(a, b)

def cn(D, g, n):
    # number of monomials of degree D in x0,x1,x2 with x2-exponent = g mod n (direct count)
    if D < 0:
        return 0
    r = g % n
    tot = 0
    b2 = r
    while b2 <= D:
        tot += D - b2 + 1
        b2 += n
    return tot

def hodge(d, n, a, chi):
    h20 = cn(d - 3 - a, chi - 1, n)
    h02 = cn(a - 3, -chi - 1, n)
    if d % n == 0:
        dim = (d - 2) * (d // n - (1 if chi % n == 0 else 0)) + (1 if (chi + a) % n == 0 else 0)
    else:
        assert (d - 1) % n == 0
        dim = (d - 1) ** 2 // n - (d - 2) * (1 if chi % n == 0 else 0)
    return h20, dim - h20 - h02, h02, dim

def group(d, n):
    if d % n == 0:
        return [(t, t % n) for t in range(1, d) if gcd(t, d) == 1]
    U1 = [t for t in range(1, d) if gcd(t, d) == 1]
    U2 = [t for t in range(1, n) if gcd(t, n) == 1] if n > 1 else [0]
    return [(t1, t2) for t1 in U1 for t2 in U2]

def in_E(d, n):
    if n == d:
        return True
    if d % n == 0 and d // n == 2:
        return True
    if d % n == 0 and d // n in (3, 4) and n % 6 == 0:
        return True
    return False

def closed_form(d, n, KR):
    if (d, n) == (6, 2):
        return 92
    if (d, n) == (6, 3):
        return 210
    if (d, n) == (6, 6):
        return 3394
    if n == d:
        return 2 + (d * d - 2 * d + 2) ** 2 + (d - 1) * (d - 2) + KR
    if n == d - 1:
        return 2 * d * d - 3 * d + 4
    if d % n == 0 and 2 * n == d:
        assert (3 * d * d - d + 2) % 2 == 0
        return (3 * d * d - d + 2) // 2
    if (d % n == 0 and 3 * n <= d) or ((d - 1) % n == 0 and 2 * n <= d - 1):
        return 3 + n * (d - 1)
    return None

def analyse(d, n, rng):
    NP = (d - 1) * n
    A = np.repeat(np.arange(1, d, dtype=np.int64), n)
    CH = np.tile(np.arange(0, n, dtype=np.int64), d - 1)
    H20 = np.zeros(NP, dtype=np.int64); H11 = np.zeros(NP, dtype=np.int64)
    H02 = np.zeros(NP, dtype=np.int64); DIM = np.zeros(NP, dtype=np.int64)
    for i in range(NP):
        h20, h11, h02, dim = hodge(d, n, int(A[i]), int(CH[i]))
        assert h11 >= 0, (d, n, A[i], CH[i])
        H20[i], H11[i], H02[i], DIM[i] = h20, h11, h02, dim
    zero = DIM == 0
    pure = (~zero) & (H20 == 0) & (H02 == 0)
    G = group(d, n)
    code = (H20 * 65536 + H11) * 65536 + H02
    code_u = code.astype(np.uint64)
    w1 = rng.integers(1, 2**62, size=len(G), dtype=np.uint64)
    w2 = rng.integers(1, 2**62, size=len(G), dtype=np.uint64)
    hs1 = np.zeros(NP, dtype=np.uint64); hs2 = np.zeros(NP, dtype=np.uint64)
    sg1 = np.zeros(NP, dtype=np.uint64)
    rep = np.arange(NP, dtype=np.int64)
    zflag = (CH % n == 0)
    for j, (t1, t2) in enumerate(G):
        ia = (t1 * A) % d
        ic = (t2 * CH) % n
        idx = (ia - 1) * n + ic
        rep = np.minimum(rep, idx)
        hs1 += w1[j] * code_u[idx]
        hs2 += w2[j] * code_u[idx]
        if d % n == 0:
            sig = ia + np.where(ic == 0, n, ic)
            sg1 += w1[j] * sig.astype(np.uint64)
    # pure orbits: orbit (canonical rep = min index) all of whose members are pure
    nonpure_member = np.zeros(NP, dtype=bool)
    np.logical_or.at(nonpure_member, rep, ~pure)
    pure_orbit = (~zero) & (~nonpure_member[rep])
    np_orbit = (~zero) & (~pure_orbit)
    NF = int(DIM[pure_orbit].sum())
    P = int(np_orbit.sum())
    # Hodge-signature groups among non-pure-orbit pieces
    sel = np.nonzero(np_orbit)[0]
    groups = {}
    for i in sel:
        groups.setdefault((int(hs1[i]), int(hs2[i])), []).append(int(i))
    fail_pairs = set()
    exact_ok = True
    for key, mem in groups.items():
        us = set(int(A[i]) for i in mem)
        if len(us) < 2:
            continue
        # exact verification of the hash group
        vecs = []
        for i in mem:
            v = []
            for (t1, t2) in G:
                ia = (t1 * int(A[i])) % d
                ic = (t2 * int(CH[i])) % n
                v.append(int(code[(ia - 1) * n + ic]))
            vecs.append(tuple(v))
        if any(v != vecs[0] for v in vecs):
            exact_ok = False
        for i in mem:
            for j in mem:
                if A[i] != A[j]:
                    fail_pairs.add((i, j))
    res = {"d": d, "n": n, "NF": NF, "P": P, "nfail": len(fail_pairs), "inE": in_E(d, n),
           "exact_ok": exact_ok}
    res["xu_ok"] = ((len(fail_pairs) > 0) == in_E(d, n))
    if d % n == 0:
        # sigma criterion (Lemmas 2.2-2.3): z = 0 and equal sigma-signature, cross-u, non-pure orbits
        sgroups = {}
        for i in sel:
            if zflag[i]:
                continue
            sgroups.setdefault(int(sg1[i]), []).append(int(i))
        spairs = set()
        for key, mem in sgroups.items():
            for i in mem:
                for j in mem:
                    if A[i] != A[j]:
                        spairs.add((i, j))
        res["sigma_agrees"] = (spairs == fail_pairs)
        z1fail = sum(1 for (i, j) in fail_pairs if zflag[i] or zflag[j])
        res["z1fail"] = z1fail
    KR = len(fail_pairs) if n == d else 0
    if n == d:
        idx = lambda a, c: (a - 1) * n + c
        swaps = set()
        for i in sel:
            a, c = int(A[i]), int(CH[i])
            if a != c:
                j = idx(c, a)
                if np_orbit[j]:
                    swaps.add((int(i), j))
        res["nswap"] = len(swaps)
        res["swaps_in_KR"] = swaps <= fail_pairs
        res["kappa"] = len(fail_pairs) - len(swaps)
        res["swapformula"] = (len(swaps) == ((d - 1) * (d - 3) + (1 if d % 2 == 0 else 0) if d >= 7 else 8))
    tot = 2 + (1 + NF) ** 2 + P + KR
    cf = closed_form(d, n, KR)
    res["dimHdg"] = tot
    res["closed"] = cf
    res["dim_ok"] = (cf == tot)
    # Theorem C range: n | d, n >= 3, m >= 2: dim != 2k and h11 >= k (consistency with ADE limit)
    if d % n == 0 and n >= 3 and d // n >= 2:
        m = d // n
        k = m * (m - 1) // 2 * n
        nz = ~zero
        res["dim_ne_2k"] = bool(np.all(DIM[nz] != 2 * k))
        res["min_h11_minus_k"] = int((H11[nz] - k).min())
    # Remark 4.6 lifts for n = d/2 (comp claim)
    if d % n == 0 and 2 * n == d and d <= 60 and fail_pairs:
        def sigsig(a, c):
            return tuple((t * a) % d + ((t * c) % d) for t in range(1, d) if gcd(t, d) == 1)
        nolift = 0
        for (i, j) in fail_pairs:
            a, c, b, c2 = int(A[i]), int(CH[i]), int(A[j]), int(CH[j])
            ok = False
            for cc in (c, c + n):
                for cc2 in (c2, c2 + n):
                    if cc % d == 0 or cc2 % d == 0:
                        continue
                    if (a + cc) % d == 0 or (b + cc2) % d == 0:
                        continue
                    if sigsig(a, cc) == sigsig(b, cc2):
                        ok = True
            if not ok:
                nolift += 1
        res["lift_missing"] = nolift
    return res, fail_pairs, (A, CH)

def main_full(dmin, dmax, out, mlist=None):
    rng = np.random.default_rng(12345)
    with open(out, "w") as f:
        for d in range(dmin, dmax + 1):
            ns = []
            for n in range(2, d + 1):
                if d % n == 0 or (d - 1) % n == 0:
                    if mlist is not None:
                        if d % n != 0 or (d // n) not in mlist:
                            continue
                    ns.append(n)
            for n in ns:
                t0 = time.time()
                res, fp, _ = analyse(d, n, rng)
                res["sec"] = round(time.time() - t0, 2)
                f.write(repr(res) + "\n")
                f.flush()
        f.write("DONE\n")

def check_pair_hodge(d, n, s, t):
    a, c = s; b, c2 = t
    assert d % n == 0
    nonpure_s = False; nonpure_t = False
    for g in range(1, d):
        if gcd(g, d) != 1:
            continue
        hs = hodge(d, n, (g * a) % d, (g * c) % n)
        ht = hodge(d, n, (g * b) % d, (g * c2) % n)
        if hs != ht:
            return "separated"
        if hs[3] == 0:
            return "zero"
        if hs[0] + hs[2] > 0:
            nonpure_s = True
        if ht[0] + ht[2] > 0:
            nonpure_t = True
    if a % d == b % d:
        return "same_u"
    if not (nonpure_s and nonpure_t):
        return "pure_orbit"
    return "FAILS_XU"

def main_fam(dmax, out):
    bad = 0; cnt = 0
    with open(out, "w") as f:
        for d in range(6, dmax + 1, 2):
            n = d // 2
            if n % 2 == 1:
                s, t = (n, (n - 1) // 2), (n - 1, (n + 1) // 2)
            elif (n // 2) % 2 == 1:
                n1 = n // 2
                s, t = (n, n1 - 1), (n - 2, n1 + 1)
            else:
                s, t = (n, n // 4), (n // 2, 3 * n // 4)
            r = check_pair_hodge(d, n, s, t)
            cnt += 1
            if r != "FAILS_XU":
                bad += 1
            f.write("A4 d=%d n=%d s=%s t=%s %s\n" % (d, n, s, t, r))
        for m in (3, 4):
            for d in range(6, dmax + 1):
                if d % m or (d // m) % 6:
                    continue
                n = d // m
                s = (d // 2, n // 6)
                t = (d // 3, (2 * n // 3) if m == 3 else (5 * n // 6))
                r = check_pair_hodge(d, n, s, t)
                cnt += 1
                if r != "FAILS_XU":
                    bad += 1
                f.write("A3 d=%d n=%d m=%d s=%s t=%s %s\n" % (d, n, m, s, t, r))
        # Remark 3.5 family at n = d (Hodge check)
        for d in range(8, dmax + 1, 2):
            s, t = (1, 1 + d // 2), (2, d // 2)
            r = check_pair_hodge(d, d, s, t)
            cnt += 1
            if r != "FAILS_XU":
                bad += 1
            f.write("R35 d=%d s=%s t=%s %s\n" % (d, s, t, r))
        for (d, s, t) in [(8, (4, 3), (6, 1)), (9, (1, 4), (2, 3)), (8, (1, 5), (2, 4))]:
            r = check_pair_hodge(d, d, s, t)
            f.write("EX d=%d s=%s t=%s %s\n" % (d, s, t, r))
        r = check_pair_hodge(8, 4, (4, 3), (6, 1))
        f.write("EX (8,4) ((4,3),(6,1)) %s\n" % r)
        f.write("TOTAL %d pairs, bad %d\nDONE\n" % (cnt, bad))

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "full":
        ml = None
        if len(sys.argv) > 5:
            ml = [int(x) for x in sys.argv[5].split(",")]
        main_full(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], ml)
    elif mode == "fam":
        main_fam(int(sys.argv[2]), sys.argv[3])
