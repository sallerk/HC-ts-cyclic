# Referee U: execute the explicit moves of the proof of Theorem 4.9, Step 3 (cases 3a, 3b, 3c),
# for EVERY ordered pair (b, b') of distinct level-0 basis monomials of every piece with h20 > 0
# (weight 2) or two-sided (weight 1), and check that the conditions of (R1)/(R2) hold exactly
# as claimed (vanishing of x^b x^m, non-wrap of x^{b'} x^m, m in S, ...).
# Usage: python -I -B refU_proofmoves.py DMIN DMAX OUTFILE
import sys, json
from multiprocessing import Pool

def in_class(d, n):
    return 2 <= n <= d and (d % n == 0 or (d - 1) % n == 0)

def make_nf(d, n):
    if d % n == 0:
        def nf(e):
            if min(e) < 0: raise ValueError
            return (1, e) if max(e) <= d - 2 else None
        return 'fermat', d - 2, nf
    def nf(e):
        if min(e) < 0: raise ValueError
        i, j, k = e; c = 1
        while True:
            if j >= d - 1 or i >= d: return None
            if i >= 1 and k >= d - 2: return None
            if i == 0 and k >= d - 1:
                c = -d * c; i = d - 1; k -= d - 1; continue
            return (c, (i, j, k))
    return 'chain', d - 1, nf

def check_pair(d, n, kind, Imax, nf, b, bp, w1, inS, nonwrap):
    i, j, k = b
    s = i + j
    D = (bp[0] - i, bp[1] - j, bp[2] - k)
    Delta = D[2]; d0 = D[0]
    sp = bp[0] + bp[1]
    if Delta % n:
        return "Delta not mult of n"
    if Delta < 0:
        aD = -Delta
        lo, hi = (d - 1 - k, d - 2 - k + aD) if kind == 'fermat' else (d - 2 - k, d - 3 - k + aD)
        m2 = next((x for x in range(lo, hi + 1) if x % n == 0), None)
        if m2 is None or not (0 <= m2 <= d - 2):
            return "3a m2"
        cap0, cap1 = Imax - bp[0], d - 2 - bp[1]
        low0 = 1 if kind == 'chain' else 0
        m0 = max(low0, d - m2 - cap1)
        if m0 > cap0 or m0 > d - m2:
            return "3a caps"
        m = (m0, d - m2 - m0, m2)
        if not inS(m): return "3a m not in S"
        if nf((i + m[0], j + m[1], k + m[2])) is not None: return "3a b*m != 0"
        if not nonwrap(tuple(bp[t] + m[t] for t in range(3))): return "3a b'*m wraps"
        return None
    # Delta >= 0
    if sp + 2 <= Imax:
        if d0 > -Delta:
            m0 = sp + 2 - bp[0]
            if m0 < 2: return "3b m0<2"
        elif d0 < 0:
            m0 = Imax - bp[0]
        else:
            return "3b impossible case reached"
        if not (0 <= m0 <= d): return "3b m0 range"
        m = (m0, d - m0, 0)
        if nf((i + m[0], j + m[1], k)) is not None: return "3b b*m != 0"
        if not nonwrap((bp[0] + m[0], bp[1] + m[1], bp[2])): return "3b b'*m wraps"
        return None
    # 3c
    if not (w1 and kind == 'fermat' and k == 0 and Delta == 0):
        return "3c reached outside weight-1 Fermat row 0"
    if d0 < 0:   # mirror x0 <-> x1
        b = (b[1], b[0], b[2]); bp = (bp[1], bp[0], bp[2]); i, j, k = b; d0 = bp[0] - i
    if j < d0: return "3c j<delta0"
    if j + 1 >= n:
        m1 = max(d - 1 - j, d - n - j - 1 + d0); m0 = d - n - m1
        if m0 < 0: return "3c m0<0"
        m = (m0, m1, n)
        if nf((i + m0, j + m1, n)) is not None: return "3c1 b*m != 0"
        if not nonwrap((bp[0] + m0, bp[1] + m1, n)): return "3c1 b'*m wraps"
        return None
    # j + 1 <= n - 1: R2 move into level 1 of row n, which is killed via b''
    if not (n <= d - 3): return "3c row n not two-level"
    bpp = (d - 3 - n, 0, n)
    if bpp[1] - d0 >= 0: return "3c b'' has partner"
    m0 = max(0, (d - n) - (d - 2 - j))
    m1 = d - n - m0
    if bp[0] + m0 > d - 2 or j + m1 > d - 2 or m1 < 0: return "3c2 caps"
    if not nonwrap((i + m0, j + m1, n)) or not nonwrap((bp[0] + m0, bp[1] + m1, n)):
        return "3c2 wraps"
    return None

def job(args):
    d, n = args
    kind, Imax, nf = make_nf(d, n)
    def inS(m):
        return min(m) >= 0 and sum(m) == d and m[2] % n == 0
    def nonwrap(e):
        r = nf(e)
        return r is not None and r[0] == 1 and r[1] == e
    def rows_ok(k):
        return 0 <= k <= (d - 2 if kind == 'fermat' else d - 3)
    fails = []; npairs = 0
    pieces = []
    for a in range(1, d):
        D0 = d - 3 - a
        for chi in range(n):
            rho = (chi - 1) % n
            L0 = [(i, D0 - k - i, k) for k in range(0, D0 + 1) if k % n == rho and rows_ok(k)
                  for i in range(0, D0 - k + 1) if D0 - k - i <= d - 2 and i <= Imax]
            if L0:
                pieces.append(((a, chi), L0, False))
    for chi in range(n):
        rho = (chi - 1) % n
        D0 = d - 3
        L0 = [(i, D0 - k - i, k) for k in range(0, D0 + 1) if k % n == rho and rows_ok(k)
              for i in range(0, D0 - k + 1) if D0 - k - i <= d - 2 and i <= Imax]
        # two-sided check: level 1 nonempty
        lvl1 = any(rows_ok(k) and k % n == rho and 0 <= 2 * d - 3 - k <= Imax + d - 2
                   for k in range(0, d)) or (kind == 'chain')
        if L0 and lvl1:
            pieces.append((("w1", chi), L0, True))
    for (lab, L0, w1) in pieces:
        for b in L0:
            for bp in L0:
                if b == bp: continue
                npairs += 1
                r = check_pair(d, n, kind, Imax, nf, b, bp, w1, inS, nonwrap)
                if r is not None:
                    fails.append((lab, b, bp, r))
    return {"d": d, "n": n, "kind": kind, "pieces": len(pieces), "pairs": npairs,
            "nfails": len(fails), "fails": fails[:5]}

def main():
    dmin, dmax, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    tasks = [(d, n) for d in range(dmax, dmin - 1, -1) for n in range(2, d + 1) if in_class(d, n)]
    tot = 0; totf = 0
    with Pool(6) as pool, open(out, "w") as f:
        for res in pool.imap_unordered(job, tasks):
            tot += res["pairs"]; totf += res["nfails"]
            if res["nfails"]:
                f.write(json.dumps(res) + "\n")
        f.write("pairs checked=%d failures=%d (d=%d..%d)\nDONE\n" % (tot, totf, dmin, dmax))

if __name__ == "__main__":
    main()
