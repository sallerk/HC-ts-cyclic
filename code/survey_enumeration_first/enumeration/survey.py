"""t24 S1 survey (v = 3): symmetric families of ternary degree-d forms F and the Hodge classes of
End_{mu_d}(H^2_prim(Y_F)), Y_F = {F + u^d = 0} in P^3 (weight-2 part), and of End(H^1(C_F)) (weight-1 part).

Family = monomial set S = all mu_n-invariant monomials, weights (0,w1,w2), n <= 2d; deduplicated up to permutation
of the variables (first (n,w) found = representative). Very general member = random coefficients mod p.
Pieces are indexed by classes kappa = [e + (1,1,1)] in C = Z^3/<S> (dual of G~ = {t in T^3: t^e = 1 for e in S},
the full diagonal symmetry group of the very general member incl. the scalars mu_d; (1,1,1) = residue form).
u-character a(kappa) = -|kappa| mod d. a != 0: piece of H_j, j = a-1, weight 2, Hodge numbers
dim(R_F)_{(q+1)d-4-j, kappa-1}, q = 0,1,2. a = 0: piece of H^1(C_F), weight 1, dims (R_F)_{d-3}, (R_F)_{2d-3}.
Hodge numbers: C-graded Koszul formula (exact, char 0), checked against R_F mod p (both primes).
Galois: kappa -> t kappa, t in (Z/|C|)^*.  Normalizer: permutations pi with pi(S) = S whose torus equation is
solvable for very general coefficients (exact lattice criterion: rows of (P_pi - I) lie in the row space of the
exponent matrix).  Criteria R1-R4 of the task (see survey report).
Usage: python -I -B survey.py d [only=NAME] [allub1=J1,J2] [max=K] [tag=TAG]
"""
import itertools
import json
import math
import os
import random
import sys
import time
from collections import Counter

import flint

sys.stdout.reconfigure(encoding='utf-8')
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, '..', 'outputs'))
sys.path.insert(0, HERE)
from ord3 import Setup  # noqa: E402  (local copy of D1_order3/code/ord3.py)

SUBSETS = [(), (0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)]
PERMS = list(itertools.permutations(range(3)))
UB1_NV_CAP = 1600


def monos3(n):
    if n < 0:
        return []
    return [(a, b, n - a - b) for a in range(n, -1, -1) for b in range(n - a, -1, -1)]


def apply(pi, e):
    return (e[pi[0]], e[pi[1]], e[pi[2]])


def pderiv(P, i, p):
    R = {}
    for e, c in P.items():
        if e[i]:
            f = list(e)
            f[i] -= 1
            f = tuple(f)
            R[f] = (R.get(f, 0) + c * e[i]) % p
    return {e: v for e, v in R.items() if v}


def rank(rows, p):
    if not rows or not rows[0]:
        return 0
    return flint.nmod_mat([list(r) for r in rows], p).rank()


class JR:
    """Jacobian ring R_F (3 variables), monomial normal forms in every degree 0..top (copied from B's famk.py)."""

    def __init__(self, F, p, d):
        self.p, self.d = p, d
        self.top = 3 * (d - 2)
        parts = [pderiv(F, i, p) for i in range(3)]
        self.mons, self.idx, self.std, self.smons, self.nfm = {}, {}, {}, {}, {}
        self.ok = True
        for n in range(self.top + 2):
            mons = monos3(n)
            idx = {e: k for k, e in enumerate(mons)}
            M = len(mons)
            gens = []
            for g in parts:
                for m in monos3(n - (d - 1)):
                    row = [0] * M
                    for e, c in g.items():
                        row[idx[(e[0] + m[0], e[1] + m[1], e[2] + m[2])]] = c
                    gens.append(row)
            if gens:
                Rm, rk = flint.nmod_mat(gens, p).rref()
                ent = [int(x) for x in Rm.entries()]
                rows = [ent[i * M:(i + 1) * M] for i in range(rk)]
            else:
                rows = []
            piv = [next(j for j, v in enumerate(r) if v) for r in rows]
            ps = set(piv)
            std = [j for j in range(M) if j not in ps]
            spos = {j: k for k, j in enumerate(std)}
            rn = len(std)
            nf = flint.nmod_mat(max(rn, 1), M, p)
            for j in std:
                nf[spos[j], j] = 1
            for r, j in zip(rows, piv):
                for s in std:
                    if r[s]:
                        nf[spos[s], j] = (-r[s]) % p
            self.mons[n], self.idx[n], self.std[n] = mons, idx, std
            self.smons[n] = [mons[j] for j in std]
            self.nfm[n] = nf
            if n == self.top + 1 and rn != 0:
                self.ok = False
            if n == self.top and rn != 1:
                self.ok = False
        self._nfcache = {}

    def dim(self, n):
        return len(self.std[n]) if 0 <= n <= self.top else 0

    def nf_mono(self, e):
        if e in self._nfcache:
            return self._nfcache[e]
        n = sum(e)
        if n > self.top or n < 0:
            out = {}
        else:
            j = self.idx[n][e]
            nf = self.nfm[n]
            out = {k: int(nf[k, j]) for k in range(self.dim(n)) if int(nf[k, j])}
        self._nfcache[e] = out
        return out

    def nf_poly(self, n, P):
        v = [0] * self.dim(n)
        for e, c in P.items():
            for k, x in self.nf_mono(e).items():
                v[k] = (v[k] + c * x) % self.p
        return v


def hnf3(vecs):
    rows = [list(v) for v in vecs if any(v)]
    basis = []
    for col in range(3):
        while True:
            nz = [r for r in rows if r[col] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda r: abs(r[col]))
            pv = nz[0]
            for r in nz[1:]:
                q = r[col] // pv[col]
                for k in range(3):
                    r[k] -= q * pv[k]
            rows = [r for r in rows if any(r)]
        nz = [r for r in rows if r[col] != 0]
        assert len(nz) == 1, "lattice not of full rank"
        pv = nz[0]
        rows = [r for r in rows if r is not pv]
        if pv[col] < 0:
            pv = [-x for x in pv]
        basis.append(pv)
    return basis


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def primes_for(L, start, k):
    out = []
    q = start - (start % L) + 1
    while len(out) < k:
        if q >= start and is_prime(q):
            out.append(q)
        q += L
    return out


# ---------------------------------------------------------------- polynomials mod p (lists, index = degree)
def ptrim(a, p):
    a = [x % p for x in a]
    while a and a[-1] == 0:
        a.pop()
    return a


def pmod(a, b, p):
    a, b = ptrim(a, p), ptrim(b, p)
    inv = pow(b[-1], p - 2, p)
    while len(a) >= len(b):
        c = a[-1] * inv % p
        sh = len(a) - len(b)
        for i in range(len(b)):
            a[sh + i] = (a[sh + i] - c * b[i]) % p
        a = ptrim(a, p)
    return a


def pgcd(a, b, p):
    a, b = ptrim(a, p), ptrim(b, p)
    while b:
        a, b = b, pmod(a, b, p)
    return a


class Ctx:
    """one random member mod p."""

    def __init__(self, S, d, p, seed, canon):
        self.S, self.d, self.p, self.canon = S, d, p, canon
        self.R = None
        for tr in range(3):
            rng = random.Random(1000003 * seed + 7919 * tr + p + 31 * d)
            F = {e: rng.randrange(1, p) for e in S}
            R = JR(F, p, d)
            if R.ok:
                self.F, self.R, self.trial = F, R, tr
                break
        self._pos = {}
        self.rng = random.Random(55 + 7 * seed + p)
        self._setup = {}

    def pos(self, n, c):
        key = (n, c)
        if key not in self._pos:
            R = self.R
            self._pos[key] = [k for k, e in enumerate(R.smons[n]) if self.canon(e) == c] if 0 <= n <= R.top else []
        return self._pos[key]

    def basis(self, n, c):
        return [self.R.smons[n][k] for k in self.pos(n, c)]

    def prod(self, m, v, n, c):
        pos = self.pos(n, c)
        if not pos:
            return []
        nf = self.R.nf_mono((m[0] + v[0], m[1] + v[1], m[2] + v[2]))
        return [nf.get(k, 0) for k in pos]

    def setup(self, j):
        if j not in self._setup:
            self._setup[j] = Setup(self.F, self.d, j, self.p)
        return self._setup[j]


# ---------------------------------------------------------------- R1 tests
def hull1(ctx, degs, c):
    P = ctx.basis(degs[0], c)
    rows = [ctx.prod(P[0], v, degs[1], c) for v in ctx.S]
    return rank(rows, ctx.p)


def hull_setup(ctx, j, c, narcs=12):
    St = ctx.setup(j)
    canon = ctx.canon
    ia = [r for r in range(St.nA) if canon(St.S[0][r]) == c]
    ib = [r for r in range(St.nB) if canon(St.S[1][St.std[r]]) == c]
    rng = random.Random(9090 + ctx.p + j)
    vecs = {0: [], 1: [], 2: [], 3: []}
    for r in ia:
        v = [0] * (len(ia) + len(ib))
        v[ia.index(r)] = 1
        vecs[0].append(v)
    for _ in range(narcs):
        L = {e: rng.randrange(1, ctx.p) for e in ctx.S}
        for r in ia:
            for k, (A, Bv) in enumerate(St.W({St.S[0][r]: 1}, L, 3), start=1):
                assert all(A[x] == 0 for x in range(St.nA) if canon(St.S[0][x]) != c), "A leak"
                vecs[k].append([A[x] for x in ia] + [Bv[x] for x in ib])
    ranks, rows = [], []
    for K in range(4):
        rows += vecs[K]
        ranks.append(rank(rows, ctx.p))
    return ranks, len(ia) + len(ib)


class SetupChi:
    """ord3's Taylor coefficients W_k (k <= 3) restricted to ONE character class c (monomial class, = kappa - 1) of
    piece j: valid whenever (R_F)_{3d-4-j, c} = 0 (h02 = 0 for this class), for any j. Same formulas as ord3.Setup
    (Griffiths reduction, constants -1, 1/2, -1/6), all matrices restricted to class-c monomials (the reduction is
    G-equivariant: m d_jF has class c iff m has class c + [e_j])."""

    def __init__(self, F, d, j, p, canon, c):
        self.p = p
        parts = [pderiv(F, i, p) for i in range(3)]
        D = [d - 4 - j + k * d for k in range(4)]
        S = [[e for e in monos3(Dk) if canon(e) == c] for Dk in D]
        idx = [{e: k for k, e in enumerate(Sk)} for Sk in S]
        E = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        self.S, self.idx = S, idx

        def cols(k):
            out = []
            for jj in range(3):
                cj = canon(tuple(c[t] + E[jj][t] for t in range(3)))
                for m in monos3(D[k] - (d - 1)):
                    if canon(m) == cj:
                        out.append((jj, m))
            return out

        def jac(k, cl):   # |S_k| x |cl|, column = m * d_jF
            M = flint.nmod_mat(max(len(S[k]), 1), max(len(cl), 1), p)
            for ci, (jj, m) in enumerate(cl):
                for e, cf in parts[jj].items():
                    M[idx[k][(e[0] + m[0], e[1] + m[1], e[2] + m[2])], ci] = cf
            return M

        def divm(k, cl):  # |S_{k-1}| x |cl|
            M = flint.nmod_mat(max(len(S[k - 1]), 1), max(len(cl), 1), p)
            for ci, (jj, m) in enumerate(cl):
                if m[jj]:
                    f = list(m)
                    f[jj] -= 1
                    M[idx[k - 1][tuple(f)], ci] = m[jj] % p
            return M
        # level 1: S_{D1,c} = Jac part (+) standard complement
        cl1 = cols(1)
        n1, nh = len(S[1]), len(cl1)
        J1 = jac(1, cl1)
        rows = [[int(J1[r, ci]) for r in range(n1)] for ci in range(nh)]
        if rows:
            Rm, rk = flint.nmod_mat(rows, p).rref()
            assert rk == nh, "Jacobian map not injective on the class"
            piv = {next(k for k in range(n1) if int(Rm[r, k])) for r in range(rk)}
        else:
            piv = set()
        self.std = [k for k in range(n1) if k not in piv]
        T = flint.nmod_mat(n1, n1, p)
        for ci in range(nh):
            for r in range(n1):
                T[r, ci] = J1[r, ci]
        for k2, s in enumerate(self.std):
            T[s, nh + k2] = 1
        Ti = T.inv()

        def rowsof(M, a, b, ncols):
            out = flint.nmod_mat(b - a, ncols, p)
            for r in range(a, b):
                for cc in range(ncols):
                    out[r - a, cc] = M[r, cc]
            return out
        self.RA = {1: divm(1, cl1) * rowsof(Ti, 0, nh, n1) if nh else flint.nmod_mat(len(S[0]), n1, p)}
        self.RB = {1: rowsof(Ti, nh, n1, n1)}
        self.nA, self.nB = len(S[0]), len(self.std)
        for k in (2, 3):
            cl = cols(k)
            nk = len(S[k])
            Mk = jac(k, cl)
            Rm, rk = Mk.rref()
            assert rk == nk, "class part of R_F nonzero in degree D%d" % k
            pv = [next(cc for cc in range(len(cl)) if int(Rm[r, cc])) for r in range(rk)]
            Sub = flint.nmod_mat(nk, nk, p)
            for r in range(nk):
                for a_, cc in enumerate(pv):
                    Sub[r, a_] = Mk[r, cc]
            Si = Sub.inv()
            RI = flint.nmod_mat(len(cl), nk, p)
            for a_, cc in enumerate(pv):
                for r in range(nk):
                    RI[cc, r] = Si[a_, r]
            Red = divm(k, cl) * RI
            self.RA[k] = self.RA[k - 1] * Red
            self.RB[k] = self.RB[k - 1] * Red

    def W(self, P, L, kmax):
        p = self.p
        out = []
        X = P
        inv2, inv6 = pow(2, p - 2, p), pow(6, p - 2, p)
        for k in range(1, kmax + 1):
            Y = {}
            for e1, c1 in X.items():
                for e2, c2 in L.items():
                    e = (e1[0] + e2[0], e1[1] + e2[1], e1[2] + e2[2])
                    Y[e] = (Y.get(e, 0) + c1 * c2) % p
            X = {e: v for e, v in Y.items() if v}
            v = flint.nmod_mat(len(self.S[k]), 1, p)
            for e, cf in X.items():
                v[self.idx[k][e], 0] = cf
            sc = {1: p - 1, 2: inv2, 3: (p - inv6) % p}[k]
            a, b = self.RA[k] * v, self.RB[k] * v
            out.append(([int(a[r, 0]) * sc % p for r in range(self.nA)],
                        [int(b[r, 0]) * sc % p for r in range(self.nB)]))
        return out


def hull_chi(ctx, j, c, narcs=12):
    St = SetupChi(ctx.F, ctx.d, j, ctx.p, ctx.canon, c)
    rng = random.Random(9191 + ctx.p + j)
    nA, nB = St.nA, St.nB
    vecs = {0: [], 1: [], 2: [], 3: []}
    for r in range(nA):
        v = [0] * (nA + nB)
        v[r] = 1
        vecs[0].append(v)
    for _ in range(narcs):
        L = {e: rng.randrange(1, ctx.p) for e in ctx.S}
        for r in range(nA):
            for k, (A, Bv) in enumerate(St.W({St.S[0][r]: 1}, L, 3), start=1):
                vecs[k].append(A + Bv)
    ranks, rows = [], []
    for K in range(4):
        rows += vecs[K]
        ranks.append(rank(rows, ctx.p))
    return ranks, nA + nB


def mat_from(rows, p):
    return flint.nmod_mat([list(r) for r in rows], p)


def nondeg2(ctx, degs, c, h1):
    """for EVERY nonzero omega in the 2-dim H^top: v -> theta_v omega maps T onto the next level (dim h1)."""
    p = ctx.p
    P = ctx.basis(degs[0], c)
    A = [mat_from([ctx.prod(Pk, v, degs[1], c) for v in ctx.S], p).transpose() for Pk in P]   # h1 x |S|
    r1, r2 = A[0].rank(), A[1].rank()
    if r1 < h1 or r2 < h1:
        return False, "M1/M2 ranks %d,%d < %d" % (r1, r2, h1)
    rng = ctx.rng
    polys = []
    for _ in range(3):
        Rm = mat_from([[rng.randrange(p) for _ in range(h1)] for _ in ctx.S], p)
        B1, B2 = A[0] * Rm, A[1] * Rm
        e1 = [int(x) for x in B1.entries()]
        e2 = [int(x) for x in B2.entries()]
        xs = list(range(h1 + 1))
        ys = []
        for s in xs:
            M = mat_from([[(s * e1[i * h1 + k] + e2[i * h1 + k]) % p for k in range(h1)] for i in range(h1)], p)
            ys.append(int(M.det()))
        V = mat_from([[pow(x, k, p) for k in range(h1 + 1)] for x in xs], p)
        cv = V.inv() * mat_from([[y] for y in ys], p)
        polys.append([int(cv[k, 0]) for k in range(h1 + 1)])
    g = polys[0]
    for q in polys[1:]:
        g = pgcd(g, q, p)
    if not g:
        return False, "determinants identically 0"
    return len(g) == 1, "gcd degree %d" % (len(g) - 1)


# ---------------------------------------------------------------- first-order (IVHS) block bound
def theta(ctx, degs, c, q, v):
    """columns: theta_v on level q -> level q+1 (list of columns)."""
    return [ctx.prod(m, v, degs[q + 1], c) for m in ctx.basis(degs[q], c)]


def ub1(ctx, ps, pt, target=0):
    """dim of {(X_q): X_{q+1} theta_v(s) = theta_v(t) X_q for all v in S} (first-order bound for Hom_Hdg(E_s,E_t)).
    Incremental over v; stops once the nullity is <= target (a known lower bound). None if too large."""
    hs, ht = ps['h'], pt['h']
    lev = len(hs)
    nv = sum(hs[q] * ht[q] for q in range(lev))
    if nv == 0:
        return 0
    neq = len(ctx.S) * sum(ht[q + 1] * hs[q] for q in range(lev - 1))
    if neq == 0:
        return nv
    if nv > UB1_NV_CAP:
        return None
    off = [0]
    for q in range(lev):
        off.append(off[-1] + hs[q] * ht[q])
    p = ctx.p
    basis = []
    rows = []
    for iv, v in enumerate(ctx.S):
        for q in range(lev - 1):
            if hs[q] == 0 or ht[q + 1] == 0:
                continue
            Ts = theta(ctx, ps['degs'], ps['c'], q, v)    # hs[q] columns of length hs[q+1]
            Tt = theta(ctx, pt['degs'], pt['c'], q, v)    # ht[q] columns of length ht[q+1]
            for a in range(ht[q + 1]):
                for b in range(hs[q]):
                    row = [0] * nv
                    for cc in range(hs[q + 1]):     # X_{q+1}[a, cc] * Ts[cc, b]
                        x = Ts[b][cc]
                        if x:
                            row[off[q + 1] + a * hs[q + 1] + cc] += x
                    for cc in range(ht[q]):         # - Tt[a, cc] * X_q[cc, b]
                        x = Tt[cc][a]
                        if x:
                            row[off[q] + cc * hs[q] + b] -= x
                    if any(row):
                        rows.append([x % p for x in row])
        if rows and (len(rows) >= nv or iv == len(ctx.S) - 1):
            Rm, rk = flint.nmod_mat(basis + rows, p).rref()
            ent = [int(x) for x in Rm.entries()]
            basis = [ent[i * nv:(i + 1) * nv] for i in range(rk)]
            rows = []
            if nv - rk <= target:
                return nv - rk
    return nv - len(basis)


def kermats(ctx, pc):
    """per level q: matrix rows (basis_q elt, coord_{q+1}) x columns v in S of theta_v."""
    out = {}
    h = pc['h']
    for q in range(len(h) - 1):
        if h[q] == 0 or h[q + 1] == 0:
            continue
        rows = []
        cols = [theta(ctx, pc['degs'], pc['c'], q, v) for v in ctx.S]
        for b in range(h[q]):
            for a in range(h[q + 1]):
                rows.append([cols[iv][b][a] for iv in range(len(ctx.S))])
        out[q] = rows
    return out


def ns1(ctx, pc):
    """first-order bound for the Hodge (1,1)-classes of a weight-2 piece: h11 - rank(H11 -> (+)_v H02)."""
    h = pc['h']
    if h[2] == 0 or h[1] == 0:
        return h[1]
    rows = []
    for Q in ctx.basis(pc['degs'][1], pc['c']):
        r = []
        for v in ctx.S:
            r += ctx.prod(Q, v, pc['degs'][2], pc['c'])
        rows.append(r)
    return h[1] - rank(rows, ctx.p)


# ---------------------------------------------------------------- enumeration
def fam_name(n, w):
    if max(w) < 10:
        return "cyc%d_%s" % (n, "".join(str(x) for x in w))
    return "cyc%d_%s" % (n, "-".join(str(x) for x in w))


def enumerate_families(d):
    allm = monos3(d)
    fams = {}
    for n in range(1, 2 * d + 1):
        for w1 in range(n):
            for w2 in range(w1, n):
                S = tuple(e for e in allm if (w1 * e[1] + w2 * e[2]) % n == 0)
                key = min(tuple(sorted(apply(pi, e) for e in S)) for pi in PERMS)
                if key not in fams:
                    fams[key] = (n, (0, w1, w2), S)
    return fams


# ---------------------------------------------------------------- one family
def analyze(n, w, S, d, opts):
    t0 = time.time()
    name = fam_name(n, w)
    res = {'name': name, 'n': n, 'w': list(w), 'd': d, 'nmon': len(S)}
    N = 3 * (d - 2)
    if flint.fmpz_mat([list(e) for e in S]).rank() < 3:
        res['status'] = 'exponents of rank < 3 (not smooth)'
        return res
    # quick necessary smoothness condition at the coordinate points
    if not all(any(e[i] >= d - 1 for e in S) for i in range(3)):
        res['status'] = 'singular at a coordinate point'
        return res
    H = hnf3(S)
    sizeC = H[0][0] * H[1][1] * H[2][2]

    def canon(x):
        x = list(x)
        for k in range(3):
            b = H[k]
            q = x[k] // b[k]
            if q:
                x = [x[i] - q * b[i] for i in range(3)]
        return tuple(x)
    res['C'] = {'hnf': H, 'order': sizeC, 'Gmax_order': sizeC // d}
    primes = primes_for(math.lcm(d, n, 8), 30000, 2)
    res['primes'] = primes
    ctxs = [Ctx(S, d, p, seed, canon) for seed, p in enumerate(primes)]
    if any(cx.R is None for cx in ctxs):
        res['status'] = 'not smooth (3 random members singular mod p)'
        return res
    meff = [rank([cx.R.nf_poly(d, {e: 1}) for e in S], cx.p) for cx in ctxs]
    res['moduli'] = meff[0]
    res['moduli_both'] = meff
    if meff[0] == 0:
        res['status'] = 'no moduli'
        return res
    res['status'] = 'ok'
    # ------------------------------------------------ Koszul Hilbert series (char 0) and check mod p
    cnt = [Counter(canon(e) for e in monos3(m)) for m in range(N + 2)]
    E = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    kcache = {}

    def kos(m, c):
        if (m, c) in kcache:
            return kcache[(m, c)]
        tot = 0
        for I in SUBSETS:
            mm = m - (d - 1) * len(I)
            if mm < 0:
                continue
            cc = canon(tuple(c[k] + sum(E[i][k] for i in I) for k in range(3)))
            tot += (-1) ** len(I) * cnt[mm].get(cc, 0)
        kcache[(m, c)] = tot
        return tot
    box = [(x0, x1, x2) for x0 in range(H[0][0]) for x1 in range(H[1][1]) for x2 in range(H[2][2])]
    mism = 0
    for cx in ctxs:
        for m in range(N + 2):
            got = Counter(canon(e) for e in cx.R.smons[m])
            for c in box:
                if got.get(c, 0) != kos(m, c):
                    mism += 1
    res['hilbert_mismatch'] = mism
    # ------------------------------------------------ pieces
    one = (1, 1, 1)
    pieces = {}
    for k in box:
        a = (-sum(k)) % d
        c = canon((k[0] - 1, k[1] - 1, k[2] - 1))
        if a:
            j = a - 1
            degs = [(q + 1) * d - 4 - j for q in range(3)]
        else:
            degs = [d - 3, 2 * d - 3]
        h = [kos(D, c) if 0 <= D <= N else 0 for D in degs]
        if sum(h):
            chi = sum(x * y for x, y in zip(w, k)) % n
            pieces[k] = {'a': a, 'c': c, 'degs': degs, 'h': h, 'chi': chi}
    lab_count = Counter((pc['a'], pc['chi']) for pc in pieces.values())
    for k, pc in pieces.items():
        pc['label'] = "(%d,%d)" % (pc['a'], pc['chi']) + ("" if lab_count[(pc['a'], pc['chi'])] == 1 else str(k))
    tot2 = sum(sum(pc['h']) for pc in pieces.values() if pc['a'])
    tot1 = sum(sum(pc['h']) for pc in pieces.values() if not pc['a'])
    assert tot2 == (d - 1) * (d * d - 3 * d + 3), (tot2, d)
    assert tot1 == (d - 1) * (d - 2), (tot1, d)
    # ------------------------------------------------ Galois, normalizer
    units = [t for t in range(1, sizeC) if math.gcd(t, sizeC) == 1]
    imgs = {}
    for k in pieces:
        imgs[k] = [canon(tuple(t * x for x in k)) for t in units]
        for i, k2 in enumerate(imgs[k]):
            assert k2 in pieces and sum(pieces[k2]['h']) == sum(pieces[k]['h'])
        hrev = pieces[imgs[k][units.index(sizeC - 1)]]['h']
        assert hrev == pieces[k]['h'][::-1], "conjugation symmetry fails"
    stab = {k: frozenset(i for i in range(len(units)) if imgs[k][i] == k) for k in pieces}
    orbit_of, orbits = {}, []
    for k in sorted(pieces):
        if k in orbit_of:
            continue
        o = sorted(set(imgs[k]))
        for k2 in o:
            orbit_of[k2] = len(orbits)
        orbits.append(o)
    Smat = [list(e) for e in S]
    Emat = [[e[i] for e in S] for i in range(3)]
    Sset = set(S)
    Sidx = {e: i for i, e in enumerate(S)}
    perms = []
    for pi in PERMS:
        if set(apply(pi, e) for e in S) != Sset:
            continue
        PmI = []
        for e in S:
            row = [0] * len(S)
            row[Sidx[apply(pi, e)]] += 1
            row[Sidx[e]] -= 1
            PmI.append(row)
        if flint.fmpz_mat(Emat + PmI).rank() == 3:
            perms.append(pi)
    res['normalizer_perms'] = [list(pi) for pi in perms]

    def pact(pi, k):
        return canon(apply(pi, k))
    for pi in perms:
        for k in pieces:
            assert pact(pi, k) in pieces and pieces[pact(pi, k)]['h'] == pieces[k]['h'], "perm breaks Hodge numbers"
    # ------------------------------------------------ orbit classification
    ocls = []
    for oi, o in enumerate(orbits):
        a0 = pieces[o[0]]['a']
        weight = 2 if a0 else 1
        hs = [pieces[k]['h'] for k in o]
        info = {'id': oi, 'weight': weight, 'size': len(o), 'emb': [[pieces[k]['label'], pieces[k]['h']] for k in o]}
        if weight == 2 and all(h[0] == 0 and h[2] == 0 for h in hs):
            info['cls'] = 'pure'
            ocls.append(info)
            continue
        if weight == 1 and all(h[0] * h[1] == 0 for h in hs):
            info['cls'] = 'CM'
            ocls.append(info)
            continue
        info['cls'] = 'residual'
        tried = []
        cands = [k for k in o if pieces[k]['h'][0] in (1, 2) and (weight == 1 or pieces[k]['h'][2] == 0)]
        cands.sort(key=lambda k: (pieces[k]['h'][0], sum(pieces[k]['h'])))
        for k in cands:
            pc = pieces[k]
            h = pc['h']
            if h[0] == 1:
                r = [hull1(cx, pc['degs'], pc['c']) for cx in ctxs]
                tr = {'emb': pc['label'], 'test': 'hull1', 'need': h[1], 'ranks': r}
                if r[0] != r[1]:
                    res.setdefault('prime_disagree_rigid', []).append(tr)
                if r[0] == h[1] and r[1] == h[1]:
                    tr['pass'] = True
                    tried.append(tr)
                    info['cls'] = 'rigid'
                    info['witness'] = tr
                    break
                if weight == 2:
                    jj = pc['a'] - 1
                    hc = [hull_chi(cx, jj, pc['c']) for cx in ctxs]
                    tr['hull_by_order'] = [x[0] for x in hc]
                    tr['dim'] = hc[0][1]
                    tr['hull_src'] = 'SetupChi'
                    if jj <= 1:
                        hh = [hull_setup(cx, jj, pc['c']) for cx in ctxs]
                        tr['hull_by_order_ord3'] = [x[0] for x in hh]
                        if [x[0] for x in hh] != [x[0] for x in hc] or hh[0][1] != hc[0][1]:
                            res.setdefault('hull_chi_vs_ord3', []).append(tr)
                    if hc[0][0] != hc[1][0]:
                        res.setdefault('prime_disagree_rigid', []).append(tr)
                    if hc[0][0][-1] == hc[0][1] and hc[1][0][-1] == hc[1][1]:
                        tr['pass'] = True
                        tr['test'] = 'hull<=3 (W_k, 12 arcs)'
                        tried.append(tr)
                        info['cls'] = 'rigid'
                        info['witness'] = tr
                        break
                tr['pass'] = False
                tried.append(tr)
            elif h[0] == 2 and h[1] >= 1:
                rr = [nondeg2(cx, pc['degs'], pc['c'], h[1]) for cx in ctxs]
                tr = {'emb': pc['label'], 'test': 'nondeg2', 'res': rr}
                if rr[0][0] != rr[1][0]:
                    res.setdefault('prime_disagree_rigid', []).append(tr)
                if rr[0][0] and rr[1][0]:
                    tr['pass'] = True
                    tried.append(tr)
                    info['cls'] = 'rigid'
                    info['witness'] = tr
                    break
                tr['pass'] = False
                tried.append(tr)
        info['tried'] = tried
        if info['cls'] == 'residual':
            if not cands:
                info['why'] = 'no embedding with h_top in {1,2}%s' % (' and h02=0' if weight == 2 else '')
            else:
                info['why'] = 'R1 tests failed'
            if weight == 2:
                info['ns1'] = min(ns1(ctxs[0], pieces[k]) for k in o)
        ocls.append(info)
    # Galois-rigid: no nontrivial Galois twist g of the orbit with h(E_tau) = h(E_{g tau}) for all tau. For a rigid
    # orbit this gives End_Hdg = K_O (g-semilinear Hodge maps have K-stable kernels, so are 0 or bijective), hence simple.
    for info, o in zip(ocls, orbits):
        if info['cls'] != 'rigid':
            continue
        s0 = o[0]
        gr = True
        for i in range(len(units)):
            if imgs[s0][i] == s0:
                continue
            if all(pieces[imgs[tau][i]]['h'] == pieces[tau]['h'] for tau in o):
                gr = False
                break
        info['galois_rigid'] = gr
    # ------------------------------------------------ automorphism span on a piece (stabilizer in the normalizer)
    span_cache, symctx = {}, {}
    res['stab_span_notes'] = []

    def stab_span(sigma):
        Hs = tuple(pi for pi in perms if pact(pi, sigma) == sigma)
        if len(Hs) <= 1:
            return 1
        if sigma in span_cache:
            return span_cache[sigma]
        if Hs not in symctx:
            # member with g_pi = pi exactly for all pi in Hs: coefficients constant on Hs-orbits of S
            orb = {}
            for e in S:
                if e in orb:
                    continue
                o = {apply(pi, e) for pi in Hs}
                for e2 in o:
                    orb[e2] = e
            cx = None
            for tr in range(3):
                rng = random.Random(4711 + tr + d)
                cf = {}
                F = {}
                for e in S:
                    r0 = orb[e]
                    if r0 not in cf:
                        cf[r0] = rng.randrange(1, primes[0])
                    F[e] = cf[r0]
                R = JR(F, primes[0], d)
                if R.ok:
                    cx = Ctx.__new__(Ctx)
                    cx.S, cx.d, cx.p, cx.canon, cx.R, cx.F = S, d, primes[0], canon, R, F
                    cx._pos, cx._setup = {}, {}
                    cx.rng = random.Random(1)
                    break
            symctx[Hs] = cx
        cx = symctx[Hs]
        if cx is None:
            res['stab_span_notes'].append('no smooth symmetric member for %s; span lb set to 1' % (Hs,))
            span_cache[sigma] = 1
            return 1
        pc = pieces[sigma]
        mats = []
        for pi in Hs:
            vec = []
            for D in pc['degs']:
                pos = cx.pos(D, pc['c'])
                for m in cx.basis(D, pc['c']):
                    nf = cx.R.nf_mono(apply(pi, m))
                    vec += [nf.get(k, 0) for k in pos]
            mats.append(vec)
        span_cache[sigma] = rank(mats, cx.p)
        return span_cache[sigma]

    # ------------------------------------------------ blocks
    vhs_cache, ub_cache = {}, {}
    allub = opts.get('allub1', set())

    def border(s, t):
        return frozenset(zip(imgs[s], imgs[t]))

    def galmin_ub1(s, t, target=0):
        B = border(s, t)
        key = (B, target)
        if key not in ub_cache:
            vals, best, arg = [], None, None
            for (x, y) in sorted(B):
                v = ub1(ctxs[0], pieces[x], pieces[y], target)
                vals.append(v)
                if v is not None:
                    if best is None or v < best:
                        best, arg = v, (x, y)
                    if best <= target:
                        break
            if best is not None and best < target:
                res.setdefault('inconsistent', []).append([pieces[s]['label'], pieces[t]['label'], best, target])
            # second prime on the minimising pair (needed for every zero / closure claim)
            v2 = ub1(ctxs[1], pieces[arg[0]], pieces[arg[1]], target) if arg is not None else None
            ub_cache[key] = (best, vals + ['p2:%s' % v2])
            if best is not None and v2 is not None and v2 != best:
                res.setdefault('prime_disagree_ub1', []).append([pieces[s]['label'], pieces[t]['label'], best, v2])
                ub_cache[key] = (max(best, v2), vals + ['p2:%s' % v2])
        return ub_cache[key]

    def vhs_compare(s, t):
        B = border(s, t)
        if B in vhs_cache:
            return vhs_cache[B]
        out = []
        differ_all = []
        for cx in ctxs:
            differ = False
            safe = False
            recs = []
            for (x, y) in sorted(B)[:6]:
                Kx, Ky = kermats(cx, pieces[x]), kermats(cx, pieces[y])
                for q in Kx:
                    r1, r2 = rank(Kx[q], cx.p), rank(Ky[q], cx.p)
                    r12 = rank(Kx[q] + Ky[q], cx.p)
                    mx1, mx2 = min(len(Kx[q]), meff[0]), min(len(Ky[q]), meff[0])
                    recs.append([pieces[x]['label'], pieces[y]['label'], q, r1, r2, r12])
                    if r12 > r1 or r12 > r2:
                        differ = True
                        if (r12 > r1 and r1 == mx1) or (r12 > r2 and r2 == mx2):
                            safe = True
            out.append({'p': cx.p, 'differ': differ, 'rank_maximal': safe, 'ranks': recs[:8]})
            differ_all.append(differ)
        vhs_cache[B] = (all(differ_all), out)
        return vhs_cache[B]

    def classify_block(s, t):
        o1, o2 = orbit_of[s], orbit_of[t]
        c1, c2 = ocls[o1]['cls'], ocls[o2]['cls']
        ps, pt = pieces[s], pieces[t]
        ds, dt = sum(ps['h']), sum(pt['h'])
        weight = 2 if ps['a'] else 1
        B = border(s, t)
        g1 = stab[s] <= stab[t]
        g2 = stab[t] <= stab[s]
        aut = [list(pi) for pi in perms if pact(pi, s) == t]
        r = {'s': ps['label'], 't': pt['label'], 'hs': ps['h'], 'ht': pt['h'], 'cls': [c1, c2]}
        if weight == 2 and 'pure' in (c1, c2):
            if c1 == c2 == 'pure':
                r.update(kind='NSxNS', lb=ds * dt, ub=ds * dt)
                return r
            other = c2 if c1 == 'pure' else c1
            if other == 'rigid':
                r.update(kind='zero', lb=0, ub=0, why='pure vs rigid (rigid orbit has no Hodge classes)')
                return r
            dp = ds if c1 == 'pure' else dt
            oo = o2 if c1 == 'pure' else o1
            r.update(kind='NSxNS-res', lb=0, ub=dp * ocls[oo].get('ns1', 0),
                     why='pure x residual: NS (x) Hdg(residual piece), algebraic')
            return r
        if weight == 1 and c1 == c2 == 'CM':
            match = all((pieces[x]['h'][0] > 0) == (pieces[y]['h'][0] > 0) for x, y in B)
            if match:
                r.update(kind='CM', lb=ds * dt, ub=ds * dt, why='CM types agree along the block orbit')
            else:
                r.update(kind='zero', lb=0, ub=0, why='CM types differ along the block orbit')
            return r
        rig = ('rigid',)
        fsrc, ftgt = {}, {}
        for x, y in sorted(B):
            fsrc.setdefault(x, []).append(y)
            ftgt.setdefault(y, []).append(x)

        def hsum(lst):
            return [sum(pieces[z]['h'][q] for z in lst) for q in range(len(ps['h']))]

        def gt(u, v):
            return any(a_ > b_ for a_, b_ in zip(u, v))
        tests = []
        if c1 in rig:
            # (i) x -> (phi(k_i x))_i is injective (its kernel is a K_O-stable rational sub-HS), so each E_x embeds
            # in the sum of its targets
            tests += [('source rigid (i): E_%s not <= sum of targets' % pieces[x]['label'], gt(pieces[x]['h'], hsum(ys)))
                      for x, ys in fsrc.items()]
        if c2 in rig:
            # (iii) the K_O'-span of im(phi) is a K_O'-stable rational sub-HS, hence all of V_O'
            tests += [('target rigid (iii): E_%s not covered by its sources' % pieces[y]['label'], gt(pieces[y]['h'], hsum(xs)))
                      for y, xs in ftgt.items()]
        bad = [m for m, b in tests if b]
        if bad:
            r.update(kind='zero', lb=0, ub=0, why='Hodge numbers: ' + bad[0])
            return r
        if c1 in rig and c2 in rig and not (g1 and g2) and ocls[o1].get('galois_rigid') and ocls[o2].get('galois_rigid'):
            r.update(kind='zero', lb=0, ub=0, why='both orbits simple (End = K, no Galois twist with equal Hodge numbers)'
                                                 ' and block orbit not bijective')
            return r
        if aut:
            sp = stab_span(s)
            if c1 in rig and c2 in rig:
                r.update(kind='aut', lb=1, ub=1, perms=aut)
                if sp != 1:
                    r['flag'] = 'rigid piece with non-scalar stabilizer span %d' % sp
                return r
            r['perms'] = aut
            r.update(kind='residual', lb=sp, why='automorphism block with non-rigid side (stabilizer span %d)' % sp)
        elif c1 in rig and c2 in rig:
            if g1 and g2:
                dif, rec = vhs_compare(s, t)
                r['vhs'] = rec
                if dif:
                    r.update(kind='zero', lb=0, ub=0, why='R3: Higgs-kernel subspaces differ (both primes)')
                    return r
                r.update(kind='residual', lb=0, why='R3: equal Hodge numbers and Higgs kernels, no automorphism')
            else:
                r.update(kind='residual', lb=0, why='stabilizers incomparable (criterion not applicable)')
        else:
            r.update(kind='residual', lb=0, why='residual orbit involved')
        u, vals = galmin_ub1(s, t, r['lb'])
        r['ub'] = u
        r['ub1_orbit'] = vals
        if u is not None and u == r['lb']:
            r['kind'] = 'aut-IVHS' if u > 0 else 'zero-IVHS'
            r['why'] += ' -> closed by first-order bound ub1 = lb = %d' % u
        return r

    groups = {}
    for k, pc in pieces.items():
        groups.setdefault(pc['a'], []).append(k)
    blocks = {}
    totals = {}
    for a in sorted(groups):
        ks = sorted(groups[a])
        bl = []
        for s in ks:
            for t in ks:
                b = classify_block(s, t)
                if a in allub and 'ub1_orbit' not in b:
                    u, vals = galmin_ub1(s, t)
                    b['ub1_check'] = u
                    b['ub1_raw'] = ub1(ctxs[0], pieces[s], pieces[t])
                elif a in allub:
                    b['ub1_check'] = b['ub']
                    b['ub1_raw'] = ub1(ctxs[0], pieces[s], pieces[t])
                bl.append(b)
        blocks[a] = bl
        lb = sum(b['lb'] for b in bl)
        ubs = [b['ub'] for b in bl]
        ub = None if any(x is None for x in ubs) else sum(ubs)
        nres = sum(1 for b in bl if b['kind'] == 'residual')
        nresR = sum(1 for b in bl if b['kind'] in ('residual', 'aut-IVHS', 'zero-IVHS'))
        totals[a] = {'lb': lb, 'ub': ub, 'n_residual_blocks': nres, 'n_residual_R': nresR, 'n_blocks': len(bl),
                     'kinds': dict(Counter(b['kind'] for b in bl))}
        if a in allub:
            totals[a]['ub1_galmin_sum'] = sum(b['ub1_check'] for b in bl) if all(
                b['ub1_check'] is not None for b in bl) else None
            totals[a]['ub1_raw_sum'] = sum(b['ub1_raw'] for b in bl) if all(
                b['ub1_raw'] is not None for b in bl) else None
    # ------------------------------------------------ summary
    Y = [o for o in ocls if o['weight'] == 2]
    Cc = [o for o in ocls if o['weight'] == 1]
    res['orbits_Y'] = dict(Counter(o['cls'] for o in Y))
    res['orbits_C'] = dict(Counter(o['cls'] for o in Cc))
    resY = sum(totals[a]['n_residual_blocks'] for a in totals if a)
    res['fully_proved_Y'] = (resY == 0)
    res['fully_proved_Y_R'] = sum(totals[a]['n_residual_R'] for a in totals if a) == 0
    res['Y_lb'] = sum(totals[a]['lb'] for a in totals if a)
    res['Y_ub'] = None if any(totals[a]['ub'] is None for a in totals if a) else sum(totals[a]['ub'] for a in totals if a)
    if 0 in totals:
        res['C_lb'], res['C_ub'] = totals[0]['lb'], totals[0]['ub']
        res['C_decided'] = totals[0]['n_residual_blocks'] == 0
    else:
        res['C_lb'] = res['C_ub'] = 0
        res['C_decided'] = True
    res['totals'] = {str(a): v for a, v in totals.items()}
    res['orbits'] = ocls
    res['pieces'] = {pc['label']: pc['h'] for pc in pieces.values()}
    res['blocks'] = {str(a): bl for a, bl in blocks.items()}
    # residual block orbits (unique)
    seen = set()
    rlist = []
    for a in sorted(blocks):
        for b in blocks[a]:
            if b['kind'] != 'residual':
                continue
            key = (min(b['s'], b['t']), max(b['s'], b['t']), a)
            if key in seen:
                continue
            seen.add(key)
            rlist.append({'a': a, 's': b['s'], 't': b['t'], 'hs': b['hs'], 'ht': b['ht'], 'cls': b['cls'],
                          'lb': b['lb'], 'ub1': b['ub'], 'why': b['why']})
    res['residual_blocks'] = rlist
    res['time_s'] = round(time.time() - t0, 1)
    return res


def fmt(res):
    if res.get('status') != 'ok':
        return "%-14s d=%d |S|=%d  %s" % (res['name'], res['d'], res['nmon'], res.get('status'))
    oy = res['orbits_Y']
    oc = res['orbits_C']
    lines = ["%-14s d=%d |S|=%d moduli=%d |C|=%d |Gmax|=%d perms=%s hilbert_mismatch=%d primes=%s time=%.0fs" % (
        res['name'], res['d'], res['nmon'], res['moduli'], res['C']['order'], res['C']['Gmax_order'],
        ["".join(map(str, p)) for p in res['normalizer_perms']], res['hilbert_mismatch'], res['primes'],
        res['time_s'])]
    lines.append("   Y orbits: rigid %d pure %d residual %d | FULLY PROVED (Y) by R1-R3: %s, with IVHS closure: %s"
                 " | Y Hodge dim: lb %s ub %s" % (
                     oy.get('rigid', 0), oy.get('pure', 0), oy.get('residual', 0),
                     'YES' if res['fully_proved_Y_R'] else 'no', 'YES' if res['fully_proved_Y'] else 'no',
                     res['Y_lb'], res['Y_ub']))
    lines.append("   C orbits: rigid %d CM %d residual %d | End_Hdg(H^1(C_F)) dim: lb %s ub %s %s" % (
        oc.get('rigid', 0), oc.get('CM', 0), oc.get('residual', 0), res['C_lb'], res['C_ub'],
        '(decided)' if res['C_decided'] else '(ub = first order, [inf])'))
    for a, t in sorted(res['totals'].items(), key=lambda x: int(x[0])):
        lines.append("   a=%s: lb %s ub %s residual %d/%d %s%s" % (
            a, t['lb'], t['ub'], t['n_residual_blocks'], t['n_blocks'], t['kinds'],
            (" ub1-check galmin %s raw %s" % (t.get('ub1_galmin_sum'), t.get('ub1_raw_sum')))
            if 'ub1_raw_sum' in t else ""))
    if res.get('inconsistent'):
        lines.append("   !! INCONSISTENT ub1 < lb: %s" % res['inconsistent'][:5])
    for k in ('prime_disagree_rigid', 'prime_disagree_ub1', 'hull_chi_vs_ord3'):
        if res.get(k):
            lines.append("   !! %s: %s" % (k, res[k][:4]))
    if res.get('stab_span_notes'):
        lines.append("   notes: %s" % res['stab_span_notes'][:3])
    for o in res['orbits']:
        if o['cls'] == 'residual':
            lines.append("   residual orbit %d (w%d): %s  why: %s%s" % (
                o['id'], o['weight'], " ".join("%s%s" % (l, h) for l, h in o['emb'][:6]) +
                (" ..." if len(o['emb']) > 6 else ""), o['why'],
                (" ns1=%s" % o.get('ns1')) if 'ns1' in o else ""))
    for b in res['residual_blocks'][:30]:
        lines.append("   RES a=%d %s%s -> %s%s cls=%s lb=%s ub1=%s : %s" % (
            b['a'], b['s'], b['hs'], b['t'], b['ht'], b['cls'], b['lb'], b['ub1'], b['why']))
    if len(res['residual_blocks']) > 30:
        lines.append("   ... %d residual block classes in total" % len(res['residual_blocks']))
    return "\n".join(lines)


def main():
    d = int(sys.argv[1])
    opts = {}
    only, mx, tag = None, None, ""
    for a in sys.argv[2:]:
        if a.startswith('only='):
            only = a[5:]
        elif a.startswith('allub1='):
            opts['allub1'] = set(int(x) for x in a[7:].split(','))
        elif a.startswith('max='):
            mx = int(a[4:])
        elif a.startswith('tag='):
            tag = a[4:]
    fams = enumerate_families(d)
    items = sorted(fams.values(), key=lambda x: (x[0], x[1]))
    if only:
        sel = []
        allm = monos3(d)
        for nm in only.split(','):
            n_ = int(nm[3:nm.index('_')])
            ws = nm[nm.index('_') + 1:]
            w_ = tuple(int(x) for x in (ws.split('-') if '-' in ws else list(ws)))
            S_ = tuple(e for e in allm if sum(x * y for x, y in zip(w_, e)) % n_ == 0)
            key = min(tuple(sorted(apply(pi, e) for e in S_)) for pi in PERMS)
            print("requested %s = family of representative %s" % (nm, fam_name(*fams[key][:2])))
            sel.append((n_, w_, S_))
        items = sel
    if mx:
        items = items[:mx]
    base = os.path.join(OUT, "survey_d%d%s" % (d, tag))
    allres = []
    print("d=%d: %d distinct monomial sets (n <= %d)" % (d, len(fams), 2 * d), flush=True)
    with open(base + ".txt", "w", encoding="utf-8", newline="\n") as fh:
        for n, w, S in items:
            res = analyze(n, w, S, d, opts)
            allres.append(res)
            txt = fmt(res)
            print(txt, flush=True)
            fh.write(txt + "\n")
            fh.flush()
            with open(base + ".json", "w", encoding="utf-8", newline="\n") as fj:
                fj.write(json.dumps(allres, default=lambda o: sorted(o) if isinstance(o, (set, frozenset)) else str(o))
                         + "\n")


if __name__ == '__main__':
    main()
