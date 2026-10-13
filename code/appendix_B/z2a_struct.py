# Z2a: structure of decomposable families Y = {x0^d + u^d + H(x1,x2) = 0}, G = Z/n acting with weights (0,1,w2)
# on (x0,x1,x2) (u fixed).  Shioda-Katsura dictionary  E_(a,chi) = L_(a0,a,a2) (x) V_chi  and the two Hodge-class
# criteria (det-line / join type, and constant CM pieces).  Conventions = Z1/V1: characters are residue
# characters sum_i w_i (m_i + 1); a = u-character = m_u + 1.
# Usage: python -I -B z2a_struct.py            (prints everything; no files written)
import sys, itertools
from math import gcd
sys.stdout.reconfigure(encoding='utf-8')


def lcm(a, b):
    return a * b // gcd(a, b)


def m2(k):
    return [(i, k - i) for i in range(k + 1)] if k >= 0 else []


class Fam:
    def __init__(self, name, d, n, w2):
        self.name, self.d, self.n, self.w = name, d, n, (1, w2)
        # weight-0 monomials of degree d in (x0,x1,x2), x0 weight 0
        self.mon = [(d - i - j, i, j) for i in range(d + 1) for j in range(d + 1 - i)
                    if (i + w2 * j) % n == 0]
        self.H = [(i, j) for (k, i, j) in self.mon if k == 0]
        self.decomp = all(k == d or k == 0 for (k, i, j) in self.mon)
        self.memo = {}

    def JH(self, k, c):
        """dim of J(H)_k, monomial weight w.m = c mod n (Koszul; dH/dx_i has weight -w_i)."""
        key = (k, c)
        if key in self.memo:
            return self.memo[key]
        d, n, w = self.d, self.n, self.w
        tot = 0
        for r in range(3):
            for I in itertools.combinations(range(2), r):
                kk = k - r * (d - 1)
                if kk < 0:
                    continue
                cc = (c + sum(w[i] for i in I)) % n
                tot += (-1) ** r * sum(1 for m in m2(kk) if (w[0] * m[0] + w[1] * m[1]) % n == cc)
        self.memo[key] = tot
        return tot

    def c_of(self, chi):
        return (chi - self.w[0] - self.w[1]) % self.n

    def V(self, chi):
        """(h10, h01) of V_chi in H^1(C_H), C_H = {w^d = H} in P^2, w of weight 0."""
        d, c = self.d, self.c_of(chi)
        h10 = sum(self.JH(d - 3 - mw, c) for mw in range(d - 1))
        h01 = sum(self.JH(2 * d - 3 - mw, c) for mw in range(d - 1))
        return (h10, h01)

    def Y(self, a, chi):
        """[h20,h11,h02] of the Y-piece (a, chi) (all a0 together)."""
        d, c = self.d, self.c_of(chi)
        out = []
        for t in (1, 2, 3):
            s = 0
            for m0 in range(d - 1):
                k = t * d - 4 - m0 - (a - 1)
                s += self.JH(k, c)
            out.append(s)
        return out

    def a0(self, a, chi):
        return (-(a + chi)) % self.d

    def L(self, a, chi):
        d = self.d
        a0 = self.a0(a, chi)
        if a0 == 0 or (a0 + a) % d == 0:
            return None
        a2 = (-(a0 + a)) % d
        return (a0, a, a2)


def pL(al, d):
    return 1 if sum(al) == d else 0


def units(N):
    return [t for t in range(1, N) if gcd(t, N) == 1]


def tmul(t, al, d):
    return tuple((t * x) % d for x in al)


def dictionary_check(F, verbose=False):
    d, n = F.d, F.n
    bad, rows = 0, []
    for a in range(1, d):
        for chi in range(n):
            hy = F.Y(a, chi)
            a0 = F.a0(a, chi)
            if a0 == 0:
                pred = [0, 0, 0]
            elif (a0 + a) % d == 0:
                pred = None  # line part (pure, algebraic)
            else:
                al = F.L(a, chi)
                h10, h01 = F.V(chi)
                pred = [h10, h01, 0] if pL(al, d) else [0, h10, h01]
            if pred is None:
                ok = (hy[0] == 0 and hy[2] == 0)
            else:
                ok = (hy == pred)
            bad += (not ok)
            if a == 1 and sum(hy) > 0:
                rows.append((a, chi, hy, F.L(a, chi), F.V(chi) if F.L(a, chi) else 'lines', ok))
    return bad, rows


def det_blocks(F):
    """Blocks (a,chi) -> (a,-chi) with V_chi of rank 2, signature (1,1): det-line criterion."""
    d, n = F.d, F.n
    N = lcm(d, n)
    U = units(N)
    res = []
    for a in range(1, d):
        for chi in range(n):
            al, be = F.L(a, chi), F.L(a, (-chi) % n)
            if al is None or be is None:
                continue
            if F.V(chi) != (1, 1):
                continue
            pt, pLt, pVt = [], [], []
            for t in U:
                ta, tb, tc = tmul(t, al, d), tmul(t, be, d), (t * chi) % n
                h10 = F.V(tc)[0]
                lp = pL(tb, d) - pL(ta, d)
                vp = 1 - h10
                pt.append(lp + vp)
                pLt.append(lp)
                pVt.append(vp)
            hod = all(x == 0 for x in pt)
            sepV = all(x == 0 for x in pVt)
            sepL = all(x == 0 for x in pLt)
            m = n // gcd(n, chi)
            res.append(dict(a=a, chi=chi, tchi=(-chi) % n, al=al, be=be, hodge=hod, sepL=sepL, sepV=sepV,
                            ordchi=m, N=lcm(d, m), p=pt, sameL=(al == be)))
    return res


def cm_pieces(F):
    """Pieces whose whole Galois orbit has definite V (constant CM); blocks between them at equal a."""
    d, n = F.d, F.n
    N = lcm(d, n)
    U = units(N)
    defin = []
    for a in range(1, d):
        for chi in range(n):
            if F.L(a, chi) is None:
                continue
            h = F.V(chi)
            if sum(h) == 0:
                continue
            if all(F.V((t * chi) % n)[0] * F.V((t * chi) % n)[1] == 0 for t in U):
                defin.append((a, chi))
    blocks = []
    for (a, c1) in defin:
        for (b, c2) in defin:
            if a != b:
                continue
            al, be = F.L(a, c1), F.L(a, c2)
            ok = True
            for t in U:
                v1, v2 = F.V((t * c1) % n), F.V((t * c2) % n)
                p1 = pL(tmul(t, al, d), d) + (1 if v1[0] > 0 else 0)
                p2 = pL(tmul(t, be, d), d) + (1 if v2[0] > 0 else 0)
                if p1 != p2:
                    ok = False
            if ok:
                blocks.append(((a, c1), (a, c2), sum(F.V(c1)) * sum(F.V(c2))))
    return defin, blocks


def report(F, show_dict=True):
    print('=' * 100)
    print(f'{F.name}: d={F.d}, G=Z/{F.n}, weights (0,1,{F.w[1]}); decomposable={F.decomp}; H monomials {F.H}')
    print('V_chi (h10,h01) in H^1(C_H):', {chi: F.V(chi) for chi in range(F.n) if sum(F.V(chi))})
    g = sum(F.V(chi)[0] for chi in range(F.n))
    print('genus check: sum h10 =', g, ' expected (d-1)(d-2)/2 =', (F.d - 1) * (F.d - 2) // 2)
    bad, rows = dictionary_check(F)
    print('dictionary check, all (a,chi):', 'OK' if bad == 0 else f'{bad} MISMATCHES')
    if show_dict:
        print('  a=1 pieces: (a,chi) Y-hodge  L=(a0,a,a2)  V=(h10,h01)')
        for r in rows:
            print('   ', r)
    db = det_blocks(F)
    print('det-line candidates (a,chi)->(a,-chi), V_chi rank 2 sig (1,1):')
    for r in db:
        print('   ', f"({r['a']},{r['chi']})->({r['a']},{r['tchi']}) L {r['al']}->{r['be']} sameL={r['sameL']} "
              f"HODGE={r['hodge']} [L-part alone {r['sepL']}, V-part alone {r['sepV']}] ord(chi)={r['ordchi']} "
              f"lcm(d,ord)={r['N']} p_t={r['p']}")
    defin, cb = cm_pieces(F)
    print('constant-CM pieces (whole Galois orbit definite):', defin)
    print('  Hodge blocks among them (s, t, dim):', cb)
    return db, defin, cb


if __name__ == '__main__':
    fams = [Fam('cyc12_017', 6, 12, 7), Fam('cyc14_018', 7, 14, 8), Fam('cyc16_019', 8, 16, 9),
            Fam('cyc18_0-1-10', 9, 18, 10), Fam('cyc20_0-1-11', 10, 20, 11)]
    for F in fams:
        report(F, show_dict=(F.d <= 7))
    # Fermat surface lines n*(1,1,1,4), degree 7 (target of the Sym^2 correspondence for cyc14_018)
    print('=' * 100)
    print('Fermat surface X_7^2, characters n*(1,1,1,4): type p = 3 - sum/7  ->',
          [(nn, tuple((nn * x) % 7 for x in (1, 1, 1, 4)), 3 - sum((nn * x) % 7 for x in (1, 1, 1, 4)) // 7)
           for nn in range(1, 7)])
    # series G = Z/2d, weights (0,1,d+1), d = 6..24 (summary only)
    print('=' * 100)
    print('SERIES x0^d + H, G = Z/2d, weights (0,1,d+1):')
    for d in range(6, 25):
        F = Fam(f'cyc{2*d}_0-1-{d+1}', d, 2 * d, d + 1)
        bad, _ = dictionary_check(F)
        db = det_blocks(F)
        hod = [r for r in db if r['hodge']]
        defin, cb = cm_pieces(F)
        nonid = [b for b in cb if b[0] != b[1]]
        print(f' d={d:2d} dict={"OK" if bad == 0 else bad} det-cands={len(db)} hodge={len(hod)} '
              f'({sum(1 for r in hod if r["sepL"] and r["sepV"])} product, '
              f'{sum(1 for r in hod if not (r["sepL"] and r["sepV"]))} Weil-type; lcm {sorted(set(r["N"] for r in hod))}) '
              f'| CM pieces {len(defin)}, CM blocks {len(cb)} (off-diag {len(nonid)}), '
              f'self-dims {sorted(set(b[2] for b in cb if b[0] == b[1]))}')
        if d <= 12:
            for r in hod:
                print(f"     ({r['a']},{r['chi']})->({r['a']},{r['tchi']}) ord(chi)={r['ordchi']} sameL={r['sameL']} "
                      f"product={r['sepL'] and r['sepV']}")
