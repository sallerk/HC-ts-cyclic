# W2: checks of the explicit (XU)-failing families and of the m in {3,4,5} classification for larger d.
# Usage: python -I -B w2fam.py DMAX_FAM DMAX_M345 OUTFILE
import sys
from math import gcd
sys.stdout.reconfigure(encoding='utf-8')

def cn(D, g, n):
    if D < 0:
        return 0
    r = g % n
    if r > D:
        return 0
    K = (D - r) // n
    return (K + 1) * (D - r + 1) - n * K * (K + 1) // 2

def dimE(d, n, a, chi):
    if d % n == 0:
        return (d - 2) * (d // n - (1 if chi % n == 0 else 0)) + (1 if (chi + a) % n == 0 else 0)
    return (d - 1) ** 2 // n - (d - 2) * (1 if chi % n == 0 else 0)

def htrip(d, n, a, chi):
    h20 = cn(d - 3 - a, chi - 1, n)
    h02 = cn(a - 3, -chi - 1, n)
    dm = dimE(d, n, a, chi)
    return (h20, dm - h20 - h02, h02)

def units(L):
    return [t for t in range(1, L) if gcd(t, L) == 1]

def sig(d, n, U, p):
    return tuple(htrip(d, n, (t * p[0]) % d, (t * p[1]) % n) for t in U)

def fails(d, n, U, s, t):
    if s[0] % d == t[0] % d:
        return False
    if dimE(d, n, *s) == 0 or dimE(d, n, *t) == 0:
        return False
    S, T = sig(d, n, U, s), sig(d, n, U, t)
    nps = any(x[0] + x[2] > 0 for x in S)
    npt = any(x[0] + x[2] > 0 for x in T)
    return nps and npt and S == T

def main():
    DF, DM, outp = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    out = open(outp, 'w', encoding='utf-8')
    def P(*x):
        s = ' '.join(str(y) for y in x)
        print(s)
        out.write(s + '\n')
    bad = 0
    # 1. n = d/2 families
    cnt = 0
    for d in range(6, DF + 1, 2):
        n = d // 2
        U = units(d)
        if n % 2 == 1:
            s, t = (n, (n - 1) // 2), (n - 1, (n + 1) // 2)
        elif n % 4 == 2:
            h = n // 2
            s, t = (n, h - 1), (n - 2, h + 1)
        else:
            s, t = (n, n // 4), (n // 2, 3 * n // 4)
        ok = fails(d, n, U, s, t)
        cnt += 1
        if not ok:
            bad += 1
            P('FAMILY_HALF_FAIL', d, n, s, t)
    P('half-families checked for d in [6,%d] even: %d pairs, bad so far %d' % (DF, cnt, bad))
    # 2. sporadic families
    cnt = 0
    for d in range(18, DF + 1, 18):
        n = d // 3
        U = units(d)
        s, t = (d // 2, n // 6), (d // 3, 2 * n // 3)
        cnt += 1
        if not fails(d, n, U, s, t):
            bad += 1
            P('FAMILY_M3_FAIL', d, n, s, t)
    for d in range(24, DF + 1, 24):
        n = d // 4
        U = units(d)
        s, t = (d // 2, n // 6), (d // 3, 5 * n // 6)
        cnt += 1
        if not fails(d, n, U, s, t):
            bad += 1
            P('FAMILY_M4_FAIL', d, n, s, t)
    P('sporadic families checked up to d=%d: %d pairs, bad so far %d' % (DF, cnt, bad))
    # 3. m in {3,4,5}: failing set nonempty iff (m in {3,4} and 6 | n)
    cnt = 0
    for d in range(6, DM + 1):
        for m in (3, 4, 5):
            if d % m:
                continue
            n = d // m
            if n < 2:
                continue
            U = units(d)
            pieces = [(a, chi) for a in range(1, d) for chi in range(n) if dimE(d, n, a, chi) > 0]
            groups = {}
            for p in pieces:
                S = sig(d, n, U, p)
                if any(x[0] + x[2] > 0 for x in S):
                    groups.setdefault(S, []).append(p)
            nf = 0
            for lst in groups.values():
                us = set(p[0] for p in lst)
                if len(us) > 1:
                    # number of unordered pairs with different u
                    k = len(lst)
                    same = {}
                    for p in lst:
                        same[p[0]] = same.get(p[0], 0) + 1
                    nf += (k * (k - 1) - sum(v * (v - 1) for v in same.values())) // 2
            pred = (m in (3, 4) and n % 6 == 0)
            cnt += 1
            if (nf > 0) != pred:
                bad += 1
                P('M345_MISMATCH', d, n, m, nf)
            if nf > 0:
                P('m345 fail', d, n, m, 'unordered failing pairs =', nf)
    P('m in {3,4,5} classes checked for d <= %d: %d pairs (d,n)' % (DM, cnt))
    P('# DONE bad=%d' % bad)
    out.close()

main()
