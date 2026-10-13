# Z2a: enumeration of decomposable families Y = {x0^d + u^d + H(x1,x2) = 0}, H very general in span(M),
# M = saturated set of exponents i of x1^i x2^(d-i).  Symmetry group S_H = full diagonal stabilizer of H
# (characters = Z^2 / L, L = lattice spanned by (i, d-i), i in M).  Reports
#  (A) det-line (join-type) blocks (a,chi) -> (a,-chi) for rank-2 eigenspaces V_chi of signature (1,1),
#      split into product type (L- and V-parts Hodge separately) and Weil type (only the product Hodge);
#  (B) non-pure pieces E = L (x) V with V constant CM (whole Galois orbit of V definite) and rank V >= 2.
# Usage: python -I -B z2a_enum.py dmin dmax
import sys, itertools
from math import gcd
sys.stdout.reconfigure(encoding='utf-8')


def lcm(a, b):
    return a * b // gcd(a, b)


def egcd(a, b):
    if b == 0:
        return (a, 1, 0) if a >= 0 else (-a, -1, 0)
    g, x, y = egcd(b, a % b)
    return (g, y, x - (a // b) * y)


class Lat:
    """HNF basis (A,B),(0,C) of a full-rank lattice in Z^2 given by generators."""

    def __init__(self, gens):
        A = 0
        for (x, y) in gens:
            A = gcd(A, x)
        det = 0
        for (p, q) in itertools.combinations(gens, 2):
            det = gcd(det, p[0] * q[1] - p[1] * q[0])
        assert A > 0 and det > 0
        self.A, self.C = A, det // A
        # vector with x = A: integer combination of generators
        cur = (0, 0)
        g = 0
        for (x, y) in gens:
            gg, s, t = egcd(g, x)
            cur = (s * cur[0] + t * x, s * cur[1] + t * y)
            g = gg
        assert cur[0] == A
        self.B = cur[1] % self.C
        self.order = self.A * self.C

    def cls(self, k1, k2):
        q = k1 // self.A
        k1 -= q * self.A
        k2 -= q * self.B
        return (k1, k2 % self.C)

    def contains(self, k1, k2):
        return self.cls(k1, k2) == (0, 0)

    def neg(self, c):
        return self.cls(-c[0], -c[1])

    def mul(self, t, c):
        return self.cls(t * c[0], t * c[1])


def analyse(d, M):
    gens = [(i, d - i) for i in M]
    L = Lat(gens)
    sat = all(((i in M) == L.contains(i, d - i)) for i in range(d + 1))
    if not sat:
        return None
    decomp = not any(L.contains(i, j) for i in range(d) for j in range(d) if 0 < i + j < d)
    h10 = {}
    for m1 in range(d - 2):
        for m2 in range(d - 2 - m1):
            c = L.cls(m1 + 1, m2 + 1)
            h10[c] = h10.get(c, 0) + 1
    chars = set(h10) | {L.neg(c) for c in h10}
    V = {c: (h10.get(c, 0), h10.get(L.neg(c), 0)) for c in chars}
    assert sum(v[0] for v in V.values()) == (d - 1) * (d - 2) // 2
    N = lcm(d, L.order)
    U = [t for t in range(1, N) if gcd(t, N) == 1]

    def s(c):
        return (c[0] + c[1]) % d

    def Lch(a, c):
        a0 = (-(a + s(c))) % d
        if a0 == 0 or (a0 + a) % d == 0:
            return None
        return (a0, a, (-(a0 + a)) % d)

    def pL(al):
        return 1 if sum(al) == d else 0

    def tm(t, al):
        return tuple((t * x) % d for x in al)

    def Vt(t, c):
        return V.get(L.mul(t, c), (0, 0))

    def order(c):
        k = 1
        while L.mul(k, c) != (0, 0):
            k += 1
        return k

    swap_ok = all((d - i) in M for i in M)
    det = []
    for c, (x, y) in V.items():
        if (x, y) != (1, 1):
            continue
        for a in range(1, d):
            al, be = Lch(a, c), Lch(a, L.neg(c))
            if al is None or be is None:
                continue
            p, pl, pv = [], [], []
            for t in U:
                lp = pL(tm(t, be)) - pL(tm(t, al))
                vp = 1 - Vt(t, c)[0]
                p.append(lp + vp)
                pl.append(lp)
                pv.append(vp)
            if all(z == 0 for z in p):
                inv = swap_ok and L.cls(c[1], c[0]) == L.neg(c)
                det.append(dict(a=a, c=c, prod=all(z == 0 for z in pl) and all(z == 0 for z in pv),
                                m=order(c), N=lcm(d, order(c)), sameL=(al == be), swapinv=inv))
    cm = []
    for c, (x, y) in V.items():
        if x + y < 2:
            continue
        if not all(Vt(t, c)[0] * Vt(t, c)[1] == 0 for t in U):
            continue
        for a in range(1, d):
            al = Lch(a, c)
            if al is None:
                continue
            types = set()
            for t in U:
                v = Vt(t, c)
                types.add(pL(tm(t, al)) + (1 if v[0] > 0 else 0))
            if types != {1}:  # not pure at every Galois embedding
                if x * y == 0 and (x > 0) == (pL(al) == 1):  # non-pure at the identity embedding
                    cm.append((a, c, x + y))
    return dict(d=d, M=M, G=(L.A, L.C), order=L.order, decomp=decomp, moduli=len(M) - 2, V=V, det=det, cm=cm,
                rank2=sum(1 for v in V.values() if v == (1, 1)))


def main():
    dmin, dmax = int(sys.argv[1]), int(sys.argv[2])
    for d in range(dmin, dmax + 1):
        print('=' * 110)
        print(f'd = {d}')
        tot = 0
        for r in range(3, d + 2):
            for M in itertools.combinations(range(d + 1), r):
                if M[0] > 1 or M[-1] < d - 1:
                    continue
                R = analyse(d, list(M))
                if R is None:
                    continue
                tot += 1
                nd = len(R['det'])
                nw = sum(1 for x in R['det'] if not x['prod'])
                # 'new' = not an identity block (m = 2: V_chi = V_-chi, det line = identity) and not an
                # automorphism block (same Fermat line and the swap maps chi to -chi)
                new = [x for x in R['det'] if x['m'] > 2 and not (x['sameL'] and x['swapinv'])]
                nnew_p = sum(1 for x in new if x['prod'])
                nnew_w = sum(1 for x in new if not x['prod'])
                bad = sorted(set(x['N'] for x in new if not x['prod'] and x['N'] > 20
                                 and any(x['N'] % q == 0 for q in range(2, x['N']))))
                if nd or R['cm']:
                    print(f"  SUMMARY d={d} M={R['M']} moduli={R['moduli']} new det blocks (all a): product {nnew_p}, "
                          f"Weil {nnew_w}; Weil lcm(d,m) outside Shioda IV range: {bad}; "
                          f"non-pure CM pieces: {len(R['cm'])}")
                if nd == 0 and not R['cm']:
                    continue
                print(f" M={R['M']} |S_H|={R['order']} HNF(A,C)={R['G']} moduli={R['moduli']} "
                      f"decomposable_family={R['decomp']} rank2(1,1)-eigenspaces={R['rank2']}")
                if nd:
                    print(f"    det-line Hodge blocks: {nd} (Weil type {nw}); "
                          f"a=1 ones: {[(x['c'], 'prod' if x['prod'] else 'Weil', 'm=%d' % x['m'], 'lcm=%d' % x['N'], 'sameL' if x['sameL'] else '', 'swap-inv' if x['swapinv'] else '') for x in R['det'] if x['a'] == 1]}")
                if R['cm']:
                    print(f"    non-pure constant-CM pieces (a, chi, rank V): {R['cm']}")
        print(f' (saturated M with >= 1 modulus and distinct roots: {tot})')


if __name__ == '__main__':
    main()
