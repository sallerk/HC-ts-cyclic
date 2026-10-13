"""t24 U1: independent checks of the statements proved in THEOREM_U.md (not imported from U2; own code).

Run:  python -I -B u1check.py MODE args
  formulas DMAX OUT          (H),(P),(F) at formula level, d = 6..DMAX: Hodge numbers c_n, zero/pure/one-sided pieces vs
                             the classification of Prop. P, pure orbits, Separation, Y/C/dim by the block rules of T'
                             vs the closed formulas, and the cross-u (Y x Y) closure count.
  conn DMIN DMAX W OUT       basis-level checks at F0 (Fermat if n|d, chain if n|d-1), W worker processes:
                             basis size, Hodge numbers from the monomial basis vs formulas, and for every non-pure
                             weight-2 piece and every two-sided weight-1 piece: the monomial graph is connected
                             (= diagonal commutant is C) and (chain) some G0-weight space of E is 1-dimensional.
  comm DMIN DMAX W OUT       full first-order commutant dimension (exact rational weighted union-find over the
                             characters of the diagonal group G of R_F0) for every non-pure / two-sided piece.
Conventions: piece E_(a,chi), a = u-character in 1..d-1, chi = (x2-exponent + 1) mod n; weight-2 levels are degrees
D0 + q d (D0 = d-3-a, q = 0,1,2); weight-1 levels are degrees d-3, 2d-3. Galois: t in (Z/lcm(d,n))^*.
"""
import sys
import math
import json
import time
from fractions import Fraction
from collections import Counter, defaultdict


# ------------------------------------------------------------------ formulas
def cnt(D, g, n):
    """c_n(D,g) = #monomials of degree D in x0,x1,x2 with x2-exponent = g mod n."""
    if D < 0:
        return 0
    r = g % n
    if r > D:
        return 0
    K = (D - r) // n
    return (K + 1) * (D - r + 1) - n * K * (K + 1) // 2


def cnt_brute(D, g, n):
    if D < 0:
        return 0
    return sum(1 for b2 in range(D + 1) if (b2 - g) % n == 0 for b1 in range(D - b2 + 1))


def pairs(dmin, dmax):
    return [(d, n) for d in range(dmin, dmax + 1) for n in range(2, d + 1) if d % n == 0 or (d - 1) % n == 0]


def hodge2(d, n, a, chi):
    h20 = cnt(d - 3 - a, chi - 1, n)
    h02 = cnt(a - 3, -chi - 1, n)
    if d % n == 0:
        dim = (d - 2) * (d // n - (1 if chi % n == 0 else 0)) + (1 if (chi + a) % n == 0 else 0)
    else:
        dim = (d - 1) ** 2 // n - (d - 2) * (1 if chi % n == 0 else 0)
    return (h20, dim - h20 - h02, h02)


def hodge1(d, n, chi):
    return (cnt(d - 3, chi - 1, n), cnt(d - 3, -chi - 1, n))


def units(L):
    return [t for t in range(1, L) if math.gcd(t, L) == 1] if L > 1 else [1]


def pure_pred(d, n):
    """Prop. P: all pure (h20 = h02 = 0, nonzero) weight-2 pieces."""
    P = set()
    if n == d:
        P |= {(a, c) for a in range(1, d) for c in range(1, d) if a + c in (d - 1, d, d + 1)}
    if n == d - 1:
        P |= {(a, 0) for a in range(1, d)}
        P |= {(a, c) for a in range(1, d) for c in range(1, d - 1) if a + c in (d - 1, d)}
    if 2 * n == d:
        P |= {(n - 1, 0), (n, 0), (n + 1, 0)}
    if 2 * n == d - 1:
        P |= {(n, 0), (n + 1, 0)}
    if (d, n) == (6, 2):
        P.add((3, 0))
    return P


def pure_orbit_pred(d, n):
    if n == d:
        return pure_pred(d, n) if d == 6 else {(a, d - a) for a in range(1, d)}
    if n == d - 1:
        return {(a, 0) for a in range(1, d)}
    if 2 * n == d:
        return {(2, 0), (3, 0), (4, 0)} if d == 6 else {(n, 0)}
    if (d, n) == (6, 2):
        return {(3, 0)}
    return set()


def Y_form(d, n):
    if (d, n) == (6, 2):
        return 73
    if (d, n) == (6, 3):
        return 69
    if (d, n) == (6, 6):
        return 681
    if n == d:
        return (d - 1) * (d * d - d - 1)
    if 2 * n == d:
        return n * (d - 1) + d * (d - 2)
    return n * (d - 1)


def C_form(d, n):
    if (d, n) == (6, 6):
        return 35
    if n == d:
        return d - 1
    if n == d - 1:
        return n - 1
    return n


def dim_form(d, n):
    special = {(6, 2): 76, (6, 3): 73, (6, 5): 30, (6, 6): 717}
    if (d, n) in special:
        return special[(d, n)]
    if n == d:
        return 1 + d * (d - 1) ** 2
    if n == d - 1:
        return d * (d - 1)
    if 2 * n == d:
        return 1 + 3 * d * d // 2 - 2 * d
    return 1 + n * d


def formulas_one(d, n):
    L = d * n // math.gcd(d, n)
    U = units(L)
    out = {'d': d, 'n': n, 'fail': []}
    H = {(a, c): hodge2(d, n, a, c) for a in range(1, d) for c in range(n)}
    nz = {s for s, h in H.items() if sum(h) > 0}
    zero = set(H) - nz
    zero_pred = {(a, 0) for a in range(1, d)} if n == d else set()
    if zero != zero_pred:
        out['fail'].append('zero pieces')
    if any(min(h) < 0 for h in H.values()):
        out['fail'].append('negative h11')
    # conjugation symmetry h20(a,chi) = h02(d-a,-chi)
    for (a, c), h in H.items():
        if h != tuple(reversed(H[(d - a, (-c) % n)])):
            out['fail'].append('conj sym %s' % ((a, c),))
            break
    pure = {s for s in nz if H[s][0] == 0 and H[s][2] == 0}
    if pure != pure_pred(d, n):
        out['fail'].append('pure pieces %s vs %s' % (sorted(pure), sorted(pure_pred(d, n))))
    orb = {}
    for s in nz:
        orb[s] = frozenset(((t * s[0]) % d, (t * s[1]) % n) for t in U)
    pure_orb_pieces = {s for s in nz if orb[s] <= pure}
    if pure_orb_pieces != pure_orbit_pred(d, n):
        out['fail'].append('pure orbits %s' % sorted(pure_orb_pieces))
    # Separation: for fixed a, (h20,h02) injective on non-pure pieces
    for a in range(1, d):
        sig = [(H[(a, c)][0], H[(a, c)][2]) for c in range(n) if (a, c) in nz and (a, c) not in pure]
        if len(sig) != len(set(sig)):
            out['fail'].append('separation a=%d' % a)
    # Y by the block rules of T'
    rho = Counter()
    for s in pure_orb_pieces:
        rho[s[0]] += sum(H[s])
    Y = sum(1 for s in nz if s not in pure_orb_pieces) + sum(v * v for v in rho.values())
    # weight 1
    H1 = {c: hodge1(d, n, c) for c in range(n)}
    nz1 = {c for c in H1 if sum(H1[c]) > 0}
    zero1_pred = {0} if n in (d, d - 1) else set()
    if set(H1) - nz1 != zero1_pred:
        out['fail'].append('w1 zero')
    one = {c for c in nz1 if 0 in H1[c]}
    one_pred = {1, d - 1} if n == d else set()
    if one != one_pred:
        out['fail'].append('w1 one-sided %s' % sorted(one))
    U1 = units(n)
    orb1 = {c: frozenset((t * c) % n for t in U1) for c in nz1}
    cm = {c for c in nz1 if orb1[c] <= one}
    if cm != ({1, 5} if (d, n) == (6, 6) else set()):
        out['fail'].append('CM orbits %s' % sorted(cm))
    sig1 = [H1[c] for c in nz1 if c not in one]
    if len(sig1) != len(set(sig1)):
        out['fail'].append('w1 separation')
    C = sum(1 for c in nz1 if c not in cm) + sum(sum(H1[c]) ** 2 for c in cm)
    out.update(Y=Y, C=C, dim=1 + Y + C, npure=len(pure), npure_orb=len(pure_orb_pieces))
    if Y != Y_form(d, n) or C != C_form(d, n) or 1 + Y + C != dim_form(d, n):
        out['fail'].append('formula Y/C/dim %d %d %d' % (Y, C, 1 + Y + C))
    # cross-u blocks (Y x Y): unclosed pairs = equal Galois-signature, different u, not both in pure orbits
    Ul = sorted(U)
    sigs = defaultdict(list)
    for s in nz:
        sigs[tuple(H[((t * s[0]) % d, (t * s[1]) % n)] for t in Ul)].append(s)
    uncl = 0
    examples = []
    for grp in sigs.values():
        for s in grp:
            for t in grp:
                if s[0] != t[0] and not (s in pure_orb_pieces and t in pure_orb_pieces):
                    uncl += 1
                    if len(examples) < 3:
                        examples.append((s, t))
    out['crossu_unclosed'] = uncl
    out['crossu_examples'] = examples
    if uncl == 0:
        rhoY = 1 + sum(sum(H[s]) for s in pure_orb_pieces)
        out['dimHdg_YxY'] = 2 + rhoY ** 2 + sum(1 for s in nz if s not in pure_orb_pieces)
    return out


def run_formulas(dmax, outpath):
    t0 = time.time()
    # c_n closed form vs brute force
    bad = 0
    for n in range(1, 13):
        for D in range(-2, 40):
            for g in range(n):
                if cnt(D, g, n) != cnt_brute(D, g, n):
                    bad += 1
    lines = ['c_n closed form vs brute force (n<=12, D<=39): %d mismatches' % bad]
    nfail = 0
    stats = Counter()
    crossu_zero = []
    for d, n in pairs(6, dmax):
        r = formulas_one(d, n)
        if r['fail']:
            nfail += 1
        cls = ('n=d' if n == d else 'n=d-1' if n == d - 1 else 'n=d/2' if 2 * n == d else
               'n=(d-1)/2' if 2 * n == d - 1 else 'n|d' if d % n == 0 else 'n|d-1')
        stats[cls] += 1
        if r['crossu_unclosed'] == 0:
            crossu_zero.append((d, n))
        if d <= 40 or r['fail']:
            lines.append('d=%3d n=%3d %-10s Y=%8d C=%4d dim=%9d pure=%3d pure-in-pure-orbits=%3d crossu_unclosed=%6d %s %s'
                         % (d, n, cls, r['Y'], r['C'], r['dim'], r['npure'], r['npure_orb'], r['crossu_unclosed'],
                            ('dimHdg(YxY)=%d' % r['dimHdg_YxY']) if 'dimHdg_YxY' in r else 'ex=%s' % r['crossu_examples'][:1],
                            'FAIL %s' % r['fail'] if r['fail'] else ''))
    lines.append('pairs per class: %s' % dict(stats))
    lines.append('pairs with all cross-u blocks closed by Galois-transported Hodge numbers: %s' % crossu_zero)
    lines.append('TOTAL pairs %d, failures %d, time %.1fs' % (sum(stats.values()), nfail, time.time() - t0))
    with open(outpath, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(lines[0])
    print('\n'.join(lines[-3:]))


# ------------------------------------------------------------------ F0 model
class Model:
    def __init__(self, d, n):
        self.d, self.n = d, n
        if d % n == 0:
            self.kind = 'fermat'
        elif (d - 1) % n == 0:
            self.kind = 'chain'
        else:
            raise ValueError
        self.S = [(i, j, d - i - j) for i in range(d + 1) for j in range(d + 1 - i) if (d - i - j) % n == 0]

    def is_basis(self, m):
        d = self.d
        i, j, k = m
        if self.kind == 'fermat':
            return max(m) <= d - 2
        if j > d - 2:
            return False
        if i == 0:
            return k <= d - 2
        return i <= d - 1 and k <= d - 3

    def nf(self, m):
        """normal form of a monomial in R_F0: (coef, basis monomial) or None."""
        d = self.d
        i, j, k = m
        if self.kind == 'fermat':
            return (1, m) if (i <= d - 2 and j <= d - 2 and k <= d - 2) else None
        if j > d - 2:
            return None
        c = 1
        if i == 0 and k >= d - 1:          # x2^(d-1) = -d x0^(d-1)
            c = -d
            i = d - 1
            k -= d - 1
        if i >= d:                         # x0^d = 0
            return None
        if i >= 1 and k >= d - 2:          # x0 x2^(d-2) = 0
            return None
        return (c, (i, j, k))

    def basis_deg(self, D):
        if D < 0:
            return []
        return [(i, j, D - i - j) for i in range(D + 1) for j in range(D + 1 - i) if self.is_basis((i, j, D - i - j))]

    def piece2(self, a, chi):
        d, n = self.d, self.n
        D0 = d - 3 - a
        return [[m for m in self.basis_deg(D0 + q * d) if (m[2] + 1 - chi) % n == 0] for q in range(3)]

    def piece1(self, chi):
        d, n = self.d, self.n
        return [[m for m in self.basis_deg(D) if (m[2] + 1 - chi) % n == 0] for D in (d - 3, 2 * d - 3)]

    def components(self, levels):
        parent = {b: b for lv in levels for b in lv}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for q in range(len(levels) - 1):
            nxt = set(levels[q + 1])
            for b in levels[q]:
                for m in self.S:
                    r = self.nf((b[0] + m[0], b[1] + m[1], b[2] + m[2]))
                    if r is not None:
                        if r[1] not in nxt:
                            raise AssertionError('left the piece')
                        ra, rb = find(b), find(r[1])
                        if ra != rb:
                            parent[ra] = rb
        return len({find(x) for x in parent})

    def unit_G0_block(self, levels):
        """chain: is some G0-weight space (level q, x1-exponent j) of E one-dimensional?"""
        c = Counter((q, b[1]) for q, lv in enumerate(levels) for b in lv)
        return any(v == 1 for v in c.values())

    # characters of the diagonal group G of R_F0 (Fermat: full torus; chain: t2^(d-1) = t0^(d-1))
    def gchar(self, m):
        if self.kind == 'fermat':
            return m
        return (m[0] + m[2], m[1], m[2] % (self.d - 1))

    def gsub(self, x, y):
        if self.kind == 'fermat':
            return (x[0] - y[0], x[1] - y[1], x[2] - y[2])
        return (x[0] - y[0], x[1] - y[1], (x[2] - y[2]) % (self.d - 1))

    def gadd(self, x, y):
        if self.kind == 'fermat':
            return (x[0] + y[0], x[1] + y[1], x[2] + y[2])
        return (x[0] + y[0], x[1] + y[1], (x[2] + y[2]) % (self.d - 1))

    def commutant_dim(self, levels):
        L = len(levels)
        lev, chmap = {}, {}
        for q, lv in enumerate(levels):
            for b in lv:
                lev[b] = q
                chmap[self.gchar(b)] = b
        if len(chmap) != len(lev):
            raise AssertionError('G-characters not distinct')
        partner = {}
        for q, lv in enumerate(levels):
            for b in lv:
                cb = self.gchar(b)
                for b2 in lv:
                    partner[(self.gsub(self.gchar(b2), cb), b)] = b2
        nodes = list(partner)
        idx = {nd: i for i, nd in enumerate(nodes)}
        N = len(nodes)
        parent = list(range(N))
        ratio = [Fraction(1)] * N
        zero = [False] * N

        def find(x):
            path = []
            while parent[x] != x:
                path.append(x)
                x = parent[x]
            root = x
            acc = Fraction(1)
            for y in reversed(path):
                acc = ratio[y] * acc
                ratio[y] = acc
                parent[y] = root
            return root

        def w(x):
            r = find(x)
            return r, (ratio[x] if x != r else Fraction(1))

        def relate(x, y, c):          # f(x) = c f(y)
            rx, wx = w(x)
            ry, wy = w(y)
            if rx == ry:
                if wx != c * wy:
                    zero[rx] = True
                return
            parent[rx] = ry
            ratio[rx] = c * wy / wx
            zero[ry] = zero[ry] or zero[rx]

        def setzero(x):
            zero[find(x)] = True
        S = self.S
        Sch = [self.gchar(m) for m in S]
        for nd, i in idx.items():
            psi, b = nd
            q = lev[b]
            b2 = partner[nd]
            if q + 1 < L:
                for m in S:
                    r1 = self.nf((b[0] + m[0], b[1] + m[1], b[2] + m[2]))
                    r2 = self.nf((b2[0] + m[0], b2[1] + m[1], b2[2] + m[2]))
                    if r1 is not None and r2 is not None:
                        j = idx.get((psi, r1[1]))
                        if j is None or partner[(psi, r1[1])] != r2[1]:
                            raise AssertionError('character bookkeeping')
                        relate(j, i, Fraction(r2[0], r1[0]))
                    elif r1 is not None:
                        j = idx.get((psi, r1[1]))
                        if j is not None:
                            setzero(j)
                    elif r2 is not None:
                        setzero(i)
            if q >= 1:   # sources b0 at level q-1 with b0*m = c*b but (psi, b0) not a node
                cb = self.gchar(b)
                for m, mc in zip(S, Sch):
                    b0 = chmap.get(self.gsub(cb, mc))
                    if b0 is not None and lev[b0] == q - 1 and (psi, b0) not in idx:
                        r = self.nf((b0[0] + m[0], b0[1] + m[1], b0[2] + m[2]))
                        if r is not None and r[1] == b:
                            setzero(i)
        return len({find(i) for i in range(N) if not zero[find(i)]})


def conn_job(dn):
    d, n = dn
    t0 = time.time()
    M = Model(d, n)
    res = {'d': d, 'n': n, 'kind': M.kind, 'fail': []}
    tot = sum(len(M.basis_deg(D)) for D in range(3 * d - 5))
    if tot != (d - 1) ** 3 or len(M.basis_deg(3 * d - 6)) != 1 or M.basis_deg(3 * d - 5):
        res['fail'].append('basis size')
    nw2 = nw1 = 0
    for a in range(1, d):
        for chi in range(n):
            lv = M.piece2(a, chi)
            h = tuple(len(x) for x in lv)
            if h != hodge2(d, n, a, chi):
                res['fail'].append('hodge %s %s' % ((a, chi), h))
            if h[0] + h[2] == 0:
                continue
            nw2 += 1
            if M.components(lv) != 1:
                res['fail'].append('disconnected w2 %s' % ((a, chi),))
            if M.kind == 'chain' and not M.unit_G0_block(lv):
                res['fail'].append('no unit G0 block w2 %s' % ((a, chi),))
    for chi in range(n):
        lv = M.piece1(chi)
        h = tuple(len(x) for x in lv)
        if h != hodge1(d, n, chi):
            res['fail'].append('hodge1 %d %s' % (chi, h))
        if h[0] * h[1] == 0:
            continue
        nw1 += 1
        if M.components(lv) != 1:
            res['fail'].append('disconnected w1 %d' % chi)
        if M.kind == 'chain' and not M.unit_G0_block(lv):
            res['fail'].append('no unit G0 block w1 %d' % chi)
    res.update(nw2=nw2, nw1=nw1, sec=round(time.time() - t0, 1))
    return res


def comm_job(dn):
    d, n = dn
    t0 = time.time()
    M = Model(d, n)
    res = {'d': d, 'n': n, 'kind': M.kind, 'fail': [], 'w2': {}, 'w1': {}}
    for a in range(1, d):
        for chi in range(n):
            lv = M.piece2(a, chi)
            if len(lv[0]) + len(lv[2]) == 0:
                continue
            c = M.commutant_dim(lv)
            res['w2']['%d,%d' % (a, chi)] = c
            if c != 1:
                res['fail'].append('w2 %s commutant %d' % ((a, chi), c))
    for chi in range(n):
        lv = M.piece1(chi)
        if len(lv[0]) * len(lv[1]) == 0:
            continue
        c = M.commutant_dim(lv)
        res['w1'][str(chi)] = c
        if c != 1:
            res['fail'].append('w1 %d commutant %d' % (chi, c))
    res['sec'] = round(time.time() - t0, 1)
    return res


def gen_ok(levels):
    """Step 1 of Thm 4.9: every basis monomial c of level q+1 is b*m (non-wrap) with b a basis monomial of level q."""
    for q in range(len(levels) - 1):
        lower = levels[q]
        for c in levels[q + 1]:
            if not any(b[0] <= c[0] and b[1] <= c[1] and b[2] <= c[2] for b in lower):
                return False
    return True


def gen_job(dn):
    d, n = dn
    t0 = time.time()
    M = Model(d, n)
    res = {'d': d, 'n': n, 'kind': M.kind, 'fail': [], 'nw2': 0, 'nw1': 0}
    for a in range(1, d):
        for chi in range(n):
            lv = M.piece2(a, chi)
            if len(lv[0]) + len(lv[2]) == 0:
                continue
            res['nw2'] += 1
            # Step 0: by duality assume h20 > 0; check generation on E or on its dual E' = (d-a, -chi)
            if len(lv[0]) > 0:
                if not gen_ok(lv):
                    res['fail'].append('w2 gen %s' % ((a, chi),))
            else:
                lv2 = M.piece2(d - a, (-chi) % n)
                if not gen_ok(lv2):
                    res['fail'].append('w2 gen dual of %s' % ((a, chi),))
    for chi in range(n):
        lv = M.piece1(chi)
        if len(lv[0]) * len(lv[1]) == 0:
            continue
        res['nw1'] += 1
        if not gen_ok(lv):
            res['fail'].append('w1 gen %d' % chi)
    res['sec'] = round(time.time() - t0, 1)
    return res


def run_pool(job, dmin, dmax, W, outpath):
    import multiprocessing as mp
    todo = sorted(pairs(dmin, dmax), key=lambda x: -x[0] * x[0] * x[0] / x[1])
    t0 = time.time()
    out = []
    with mp.Pool(W) as pool:
        for r in pool.imap_unordered(job, todo):
            out.append(r)
            with open(outpath, 'a') as f:
                f.write(json.dumps(r) + '\n')
    nf = sum(1 for r in out if r['fail'])
    summ = 'DONE %s d=%d..%d pairs=%d failures=%d time=%.1fs' % (job.__name__, dmin, dmax, len(out), nf, time.time() - t0)
    with open(outpath, 'a') as f:
        f.write(summ + '\n')
    print(summ)


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'formulas':
        run_formulas(int(sys.argv[2]), sys.argv[3])
    elif mode == 'conn':
        run_pool(conn_job, int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
    elif mode == 'comm':
        run_pool(comm_job, int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
    elif mode == 'gen':
        run_pool(gen_job, int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
    elif mode == 'xu':
        # list the pairs d = 6..DMAX NOT in the 'closed' list of a formulas output file
        import ast
        dmax, path = int(sys.argv[2]), sys.argv[3]
        txt = open(path).read()
        key = 'closed by Galois-transported Hodge numbers: '
        closed = set(ast.literal_eval(txt.split(key)[1].split('\n')[0]))
        un = [p for p in pairs(6, dmax) if p not in closed]
        cls = Counter('n=d' if n == d else 'n=d/2' if 2 * n == d else 'other' for d, n in un)
        print('unclosed pairs:', len(un), dict(cls))
        print('unclosed, n != d and n != d/2:', [p for p in un if p[1] != p[0] and 2 * p[1] != p[0]])
        print('n=d or n=d/2 pairs that ARE closed:', sorted(p for p in closed if p[1] == p[0] or 2 * p[1] == p[0]))
    else:
        raise SystemExit('unknown mode')
