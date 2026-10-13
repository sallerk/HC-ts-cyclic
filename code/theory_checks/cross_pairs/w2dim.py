# W2: check of the closed forms of THEOREM_YY 5.1 (without the KR term), and of Lemma 4.3 (dim != 2k).
# Usage: python -I -B w2dim.py DMAX OUTFILE
import sys
from math import gcd, comb
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

def pure(d, n, a, chi):
    return dimE(d, n, a, chi) > 0 and cn(d - 3 - a, chi - 1, n) == 0 and cn(a - 3, -chi - 1, n) == 0

def closed(d, n):
    if (d, n) == (6, 2):
        return 92
    if (d, n) == (6, 3):
        return 210
    if (d, n) == (6, 6):
        return 2 + 58 ** 2 + 12
    if n == d:
        return 2 + (d * d - 2 * d + 2) ** 2 + (d - 1) * (d - 2)
    if n == d - 1:
        return 2 * d * d - 3 * d + 4
    if 2 * n == d:
        return (3 * d * d - d + 2) // 2
    return 3 + n * (d - 1)

def main():
    D, outp = int(sys.argv[1]), sys.argv[2]
    out = open(outp, 'w', encoding='utf-8')
    bad = 0
    cnt = 0
    for d in range(6, D + 1):
        for n in range(2, d + 1):
            if not (d % n == 0 or (d - 1) % n == 0):
                continue
            L = d * n // gcd(d, n)
            U = [t for t in range(1, L) if gcd(t, L) == 1]
            Nf = 0
            P = 0
            for a in range(1, d):
                for chi in range(n):
                    dm = dimE(d, n, a, chi)
                    if dm == 0:
                        continue
                    if all(pure(d, n, (t * a) % d, (t * chi) % n) for t in U):
                        Nf += dm
                    else:
                        P += 1
            val = 2 + (1 + Nf) ** 2 + P
            cnt += 1
            if val != closed(d, n):
                bad += 1
                out.write('MISMATCH %d %d %d %d\n' % (d, n, val, closed(d, n)))
            # Lemma 4.3: n | d, 3 <= n <= d/2: dim != 2k for every nonzero piece
            if d % n == 0 and n >= 3 and 2 * n <= d:
                m = d // n
                k = comb(m, 2) * n
                for a in range(1, d):
                    for chi in range(n):
                        dm = dimE(d, n, a, chi)
                        if dm > 0 and dm == 2 * k:
                            bad += 1
                            out.write('L43 %d %d %d %d\n' % (d, n, a, chi))
    out.write('checked %d pairs (d,n), d <= %d; bad=%d\n' % (cnt, D, bad))
    print('checked', cnt, 'bad', bad)
    out.close()

main()
