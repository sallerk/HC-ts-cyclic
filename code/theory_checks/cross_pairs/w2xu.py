# W2: (XU) classification checks.  Own code (conventions of THEOREM_U Prop. 2.3).
# Usage: python -I -B w2xu.py DMAX OUTFILE
# For every (d,n) in the class (d >= 6, 2 <= n <= d, n | d or n | d-1), d <= DMAX:
#  - Hodge triple h(a,chi) = (c_n(d-3-a,chi-1), dim - h20 - h02, c_n(a-3,-chi-1)), dim from Prop. 2.3;
#  - signature of a piece s = tuple over t in (Z/L)^* of h(t s); (XU) fails for {s,t'} iff equal signatures
#    (both in non-pure orbits, u(s) != u(t'));
#  - compares with the predicted classification of THEOREM_YY Thm A, and with the "sigma criterion"
#    (Lemma A.2 / A.4) computed independently from fractional parts;
#  - for n = d: checks failure <=> sigma_g(s) = sigma_g(t) for all g, counts KR(d) (ordered), swaps;
#  - for n = d/2: checks whether failing pairs have a lift in KR(d); tests candidate universal families.
import sys
from fractions import Fraction as Fr
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

def Nstar(x, n):
    # closed form of THEOREM_YY Lemma A.1
    if x <= 0:
        return 0
    e = x % n
    if e == 0:
        e = n
    v = x * x + n * x + e * (n - e)
    assert v % (2 * n) == 0
    return v // (2 * n)

def units(L):
    return [t for t in range(1, L) if gcd(t, L) == 1]

def inclass(d, n):
    return 2 <= n <= d and (d % n == 0 or (d - 1) % n == 0)

def frac(x):
    return x - (x.numerator // x.denominator)

def sigma_crit(d, n, s, t, U):
    # Lemma A.4 criterion (n | d, m = d/n >= 2, z = 0 for both): returns True iff predicted failure
    m = d // n
    a, chi = s
    b, chi2 = t
    if chi % n == 0 or chi2 % n == 0 or a == b:
        return False
    if (chi2 - chi - (a - b)) % n != 0:
        return False
    th = Fr(a - b, d)
    if frac(m * th) == 0:
        return False
    al = Fr(a, d)
    ga = Fr(chi, n)
    for u in U:
        ph = frac(u * th)
        if Fr(1, m) <= ph < 1 - Fr(1, m):
            return False
        if ph < Fr(1, m):
            if not (frac(u * al) > ph and frac(u * ga) + m * ph < 1):
                return False
    return True

def main():
    DMAX = int(sys.argv[1])
    out = open(sys.argv[2], 'w', encoding='utf-8')
    def P(*x):
        s = ' '.join(str(y) for y in x)
        print(s)
        out.write(s + '\n')
    bad = 0
    P('# w2xu.py DMAX=%d' % DMAX)
    P('# d n class nfail_unordered predicted_nonempty sigma_crit_agrees extra')
    kr_counts = {}
    for d in range(6, DMAX + 1):
        for n in range(2, d + 1):
            if not inclass(d, n):
                continue
            L = d * n // gcd(d, n)
            U = units(L)
            pieces = [(a, chi) for a in range(1, d) for chi in range(n) if dimE(d, n, a, chi) > 0]
            # Nstar check on all pieces (z=0 / z=1 sigma formula)
            for (a, chi) in pieces:
                r = chi % n if chi % n else n
                rt = (-chi) % n if (-chi) % n else n
                h = htrip(d, n, a, chi)
                assert h[0] == Nstar(d - a - 1 - r, n), (d, n, a, chi)
                assert h[2] == Nstar(a - 1 - rt, n), (d, n, a, chi)
                assert h[1] >= 0
            sig = {}
            nonpure_orbit = {}
            for (a, chi) in pieces:
                v = tuple(htrip(d, n, (t * a) % d, (t * chi) % n) for t in U)
                sig[(a, chi)] = v
                nonpure_orbit[(a, chi)] = any(x[0] + x[2] > 0 for x in v)
            groups = {}
            for p in pieces:
                if nonpure_orbit[p]:
                    groups.setdefault(sig[p], []).append(p)
            fails = []
            for g, lst in groups.items():
                for i in range(len(lst)):
                    for j in range(i + 1, len(lst)):
                        if lst[i][0] != lst[j][0]:
                            fails.append((lst[i], lst[j]))
            m = d // n if d % n == 0 else None
            if d % n == 0:
                cls = 'n|d m=%d' % m
            else:
                cls = 'n|d-1'
            pred = (n == d) or (2 * n == d) or (m in (3, 4) and n % 6 == 0)
            ok = (len(fails) > 0) == pred
            agree = '-'
            extra = ''
            if m is not None and m >= 2:
                # independent sigma criterion over ordered pairs
                fs = set()
                for (s, t) in fails:
                    fs.add((s, t)); fs.add((t, s))
                cs = set()
                bys = {}
                for p in pieces:
                    if nonpure_orbit[p] and p[1] % n != 0:
                        bys.setdefault(p[0] + (p[1] % n), []).append(p)
                for lst in bys.values():
                    for s in lst:
                        for t in lst:
                            if s != t and sigma_crit(d, n, s, t, U):
                                cs.add((s, t))
                agree = (fs == cs)
                if not agree:
                    bad += 1
                    extra += ' MISMATCH fs-cs=%s cs-fs=%s' % (sorted(fs - cs)[:4], sorted(cs - fs)[:4])
            if n == d:
                # failure <=> sigma_g equal for all g ; KR(d) ordered count ; swaps
                def sg(p):
                    return tuple(((t * p[0]) % d) + ((t * p[1]) % d) for t in U)
                fs = set()
                for (s, t) in fails:
                    fs.add((s, t)); fs.add((t, s))
                byS = {}
                for p in pieces:
                    if nonpure_orbit[p]:
                        byS.setdefault(sg(p), []).append(p)
                kr = set()
                for g, lst in byS.items():
                    for s in lst:
                        for t in lst:
                            if s[0] != t[0]:
                                kr.add((s, t))
                sw = set((s, (s[1], s[0])) for s in pieces if nonpure_orbit[s] and s[0] != s[1])
                assert sw <= kr
                if fs != kr:
                    bad += 1
                    extra += ' KR_MISMATCH'
                nsw = (d - 1) * (d - 3) + (1 if d % 2 == 0 else 0)
                if len(sw) != nsw:
                    bad += 1
                    extra += ' SWAPCOUNT %d %d' % (len(sw), nsw)
                kr_counts[d] = (len(kr), len(sw))
                extra += ' KR_ordered=%d swap_ordered=%d nonswap=%d' % (len(kr), len(sw), len(kr) - len(sw))
                if gcd(d, 6) == 1 and len(kr) != len(sw):
                    extra += ' NONSWAP_WITH_gcd6=1'
            if 2 * n == d:
                # lifts in KR(d)
                Ud = units(d)
                def sgd(p):
                    return tuple(((t * p[0]) % d) + ((t * p[1]) % d) for t in Ud)
                def inKR(s, t):
                    return s[0] != t[0] and sgd(s) == sgd(t) and (s[0] + s[1]) % d != 0
                nolift = 0
                bothlift = 0
                for (s, t) in fails:
                    L1 = [(s[0], s[1]), (s[0], s[1] + n)]
                    L2 = [(t[0], t[1]), (t[0], t[1] + n)]
                    hits = [(x, y) for x in L1 for y in L2 if inKR(x, y)]
                    if not hits:
                        nolift += 1
                    # perfect matching of both lifts
                    if (inKR(L1[0], L2[0]) and inKR(L1[1], L2[1])) or (inKR(L1[0], L2[1]) and inKR(L1[1], L2[0])):
                        bothlift += 1
                extra += ' half: nolift=%d bothlift=%d' % (nolift, bothlift)
            if not ok:
                bad += 1
                extra += ' PRED_MISMATCH'
            P(d, n, cls, len(fails), pred, agree, extra)
    P('# DONE bad=%d' % bad)
    out.close()

main()
