# V1 independent implementation: Jacobian rings mod p, Griffiths reduction, flat Taylor expansion,
# hull ranks, Higgs fields, L4' splitting test, Higgs kernels, first-order block bounds.
# Written from scratch from S2_lemmas.md (spec); no code of S1/S2/t23 consulted.
import flint
import itertools
from fractions import Fraction
from math import gcd


def monos(deg):
    if deg < 0:
        return []
    return [(a, b, deg - a - b) for a in range(deg, -1, -1) for b in range(deg - a, -1, -1)]


def madd(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def nm(rows, ncols, p):
    nr = len(rows)
    if nr == 0 or ncols == 0:
        raise ValueError("empty matrix")
    return flint.nmod_mat(nr, ncols, [x for r in rows for x in r], p)


def tl(M):
    nr, nc = M.nrows(), M.ncols()
    e = [int(x) for x in M.entries()]
    return [e[i * nc:(i + 1) * nc] for i in range(nr)]


def rank_rows(rows, ncols, p):
    if not rows or ncols == 0:
        return 0
    return nm(rows, ncols, p).rank()


def left_kernel(rows, ncols, p):
    """basis (rref) of {c : sum_i c_i rows[i] = 0}; rows: m vectors of length ncols"""
    m = len(rows)
    if m == 0:
        return []
    if ncols == 0:
        return [tuple(1 if i == j else 0 for i in range(m)) for j in range(m)]
    MT = nm([[rows[i][j] for i in range(m)] for j in range(ncols)], m, p)
    X, nl = MT.nullspace()
    if nl == 0:
        return []
    Xl = tl(X)
    vecs = [[Xl[i][j] for i in range(m)] for j in range(nl)]
    R, rk = nm(vecs, m, p).rref()
    return [tuple(r) for r in tl(R)[:rk]]


def koszul_dim(e, c, d, w, n):
    tot = 0
    for r in range(4):
        for I in itertools.combinations(range(3), r):
            ee = e - r * (d - 1)
            cc = (c + sum(w[i] for i in I)) % n
            cnt = sum(1 for m in monos(ee) if (w[0] * m[0] + w[1] * m[1] + w[2] * m[2]) % n == cc)
            tot += (-1) ** r * cnt
    return tot


class Ring:
    """R_F = S/(dF/dx0, dF/dx1, dF/dx2) mod p, graded by degree and Gamma-weight (w.m mod n)."""

    def __init__(self, d, w, n, p, F):
        self.d, self.w, self.n, self.p = d, tuple(w), n, p
        self.F = {m: c % p for m, c in F.items() if c % p}
        self.dF = []
        for i in range(3):
            g = {}
            for m, c in self.F.items():
                if m[i]:
                    mm = list(m)
                    mm[i] -= 1
                    mm = tuple(mm)
                    g[mm] = (g.get(mm, 0) + c * m[i]) % p
            self.dF.append({k: v for k, v in g.items() if v})
        self._cols, self._mac, self._red = {}, {}, {}

    def wt(self, m):
        return (self.w[0] * m[0] + self.w[1] * m[1] + self.w[2] * m[2]) % self.n

    def cols(self, e, c):
        key = (e, c)
        if key not in self._cols:
            L = [m for m in monos(e) if self.wt(m) == c]
            self._cols[key] = (L, {m: k for k, m in enumerate(L)})
        return self._cols[key]

    def mac(self, e, c):
        """Macaulay data in degree e, weight c: rref of [M | I], M rows = x^m dF_i."""
        key = (e, c)
        if key in self._mac:
            return self._mac[key]
        p, d, n = self.p, self.d, self.n
        cl, ci = self.cols(e, c)
        rows = [(i, m) for i in range(3) for m in monos(e - d + 1) if (self.wt(m) - self.w[i]) % n == c]
        nc, nr = len(cl), len(rows)
        if nc == 0 or nr == 0:
            data = dict(piv=[], nonpiv=list(range(nc)), rows=rows, Rrows=[], Trows=[])
        else:
            W = nc + nr
            flat = [0] * (nr * W)
            for r, (i, m) in enumerate(rows):
                base = r * W
                for t, cf in self.dF[i].items():
                    flat[base + ci[madd(m, t)]] = cf
                flat[base + nc + r] = 1
            A = flint.nmod_mat(nr, W, flat, p)
            Rr, rk = A.rref()
            ent = [int(x) for x in Rr.entries()]
            piv, Rrows, Trows = [], [], []
            for r in range(nr):
                row = ent[r * W:(r + 1) * W]
                j = next((k for k in range(W) if row[k]), None)
                if j is None or j >= nc:
                    break
                piv.append(j)
                Rrows.append(row[:nc])
                Trows.append(row[nc:])
            ps = set(piv)
            data = dict(piv=piv, nonpiv=[k for k in range(nc) if k not in ps], rows=rows, Rrows=Rrows, Trows=Trows)
        self._mac[key] = data
        return data

    def dimR(self, e, c):
        if e < 0:
            return 0
        if e < self.d - 1:
            return len(self.cols(e, c)[0])
        return len(self.mac(e, c)['nonpiv'])

    def std(self, e, c):
        if e < 0:
            return []
        cl = self.cols(e, c)[0]
        if e < self.d - 1:
            return list(cl)
        return [cl[k] for k in self.mac(e, c)['nonpiv']]

    def check_hilbert(self, emax=None):
        d = self.d
        if emax is None:
            emax = 3 * d - 5
        bad = []
        for e in range(0, emax + 1):
            for c in range(self.n):
                ex, got = koszul_dim(e, c, d, self.w, self.n), self.dimR(e, c)
                if ex != got:
                    bad.append((e, c, ex, got))
        return bad

    def red(self, shift, q, c):
        """Griffiths reduction of eta(x^m, q) := Res(u^j x^m Omega / G^q) (shift = 4+j) or
        Res(x^m Omega_2 / F^q) (shift = 3), deg x^m = q d - shift, to the basis
        {eta(M, q'), q' <= q, M standard of degree q' d - shift}.
        Returns (rows, dims): rows[i] = coordinates of eta(cols[i], q), layout level 1 | level 2 | ... | level q.
        Rule: eta(sum g_i F_i, q) = 1/(q-1) eta(sum d_i g_i, q-1)."""
        key = (shift, q, c)
        if key in self._red:
            return self._red[key]
        d, p = self.d, self.p
        e = q * d - shift
        dims = [self.dimR(qq * d - shift, c) for qq in range(1, q + 1)]
        nlow, ntop = sum(dims[:-1]), dims[-1]
        if e < 0:
            res = ([], dims)
        else:
            cl, ci = self.cols(e, c)
            nc = len(cl)
            if e < d - 1:
                rows = [[0] * nlow + [1 if k == r else 0 for k in range(nc)] for r in range(nc)]
            else:
                mc = self.mac(e, c)
                npos = {k: idx for idx, k in enumerate(mc['nonpiv'])}
                top = [None] * nc
                for k in mc['nonpiv']:
                    v = [0] * ntop
                    v[npos[k]] = 1
                    top[k] = v
                for r, k in enumerate(mc['piv']):
                    Rr = mc['Rrows'][r]
                    top[k] = [(-Rr[kk]) % p for kk in mc['nonpiv']]
                low = [[0] * nlow for _ in range(nc)]
                if nlow > 0 and mc['piv'] and q >= 2:
                    cl2, ci2 = self.cols(e - d, c)
                    if cl2:
                        lowrows, _ = self.red(shift, q - 1, c)
                        nr = len(mc['rows'])
                        D = [[0] * len(cl2) for _ in range(nr)]
                        for r, (i, m) in enumerate(mc['rows']):
                            if m[i]:
                                mm = list(m)
                                mm[i] -= 1
                                D[r][ci2[tuple(mm)]] = m[i] % p
                        inv = pow(q - 1, p - 2, p)
                        prod = nm(mc['Trows'], nr, p) * nm(D, len(cl2), p) * nm(lowrows, nlow, p)
                        pl = tl(prod)
                        for r, k in enumerate(mc['piv']):
                            low[k] = [(x * inv) % p for x in pl[r]]
                rows = [low[k] + top[k] for k in range(nc)]
            res = (rows, dims)
        self._red[key] = res
        return res

    def reduce_poly(self, shift, q, c, poly):
        rows, dims = self.red(shift, q, c)
        ci = self.cols(q * self.d - shift, c)[1]
        out = [0] * sum(dims)
        p = self.p
        for m, cf in poly.items():
            cf %= p
            if not cf:
                continue
            r = rows[ci[m]]
            for k in range(len(out)):
                if r[k]:
                    out[k] = (out[k] + cf * r[k]) % p
        return out, dims

    def hull(self, shift, c, S, Tk, qmax, kmax=3):
        d = self.d
        dims = [self.dimR(qq * d - shift, c) for qq in range(1, qmax + 1)]
        dimE = sum(dims)
        Pl = self.cols(d - shift, c)[0]
        allrows, ranks = [], []
        for k in range(kmax + 1):
            q = k + 1
            rows_q, dq = self.red(shift, q, c)
            ci = self.cols(q * d - shift, c)[1]
            width = sum(dq)
            for P in Pl:
                for m in Tk[k]:
                    r = rows_q[ci[madd(P, m)]]
                    if width < dimE:
                        r = r + [0] * (dimE - width)
                    elif width > dimE:
                        assert not any(r[dimE:])
                        r = r[:dimE]
                    allrows.append(r)
            rk = rank_rows(allrows, dimE, self.p)
            ranks.append(rk)
            if rk == dimE:
                ranks += [dimE] * (kmax - k)
                break
        return ranks, dimE

    def theta(self, shift, c, S, qmax):
        """thetas[l][vi] = matrix (h_{l+1} x h_l) of theta_v : Gr level l -> level l+1 (up to the factor -l)."""
        d = self.d
        out = []
        for l in range(1, qmax):
            src = self.std(l * d - shift, c)
            rows_n, dn = self.red(shift, l + 1, c)
            ci = self.cols((l + 1) * d - shift, c)[1]
            off, ht = sum(dn[:-1]), dn[-1]
            mats = []
            for v in S:
                M = [[0] * len(src) for _ in range(ht)]
                for kc, P in enumerate(src):
                    r = rows_n[ci[madd(P, v)]]
                    for i in range(ht):
                        M[i][kc] = r[off + i]
                mats.append(M)
            out.append(mats)
        return out


def interp(xs, ys, p):
    x = flint.nmod_poly([0, 1], p)
    res = flint.nmod_poly([], p)
    for i, xi in enumerate(xs):
        num = flint.nmod_poly([1], p)
        den = 1
        for j, xj in enumerate(xs):
            if j != i:
                num = num * (x - xj)
                den = den * (xi - xj) % p
        res = res + num * (int(ys[i]) * pow(den, p - 2, p) % p)
    return res


def minrank_cert(M1, M2, k, p, rng, tries=4):
    """True => rank(s M1 + t M2) >= k for every (s:t) in P^1(F_p-bar) (certificate: gcd of compressed
    k x k determinants is constant and the point at infinity is excluded). False = no certificate."""
    h = len(M1)
    m = len(M1[0]) if h else 0
    if k <= 0:
        return True
    if k > min(h, m):
        return False
    g = None
    tops = []
    for _ in range(tries):
        P = nm([[rng.randrange(p) for _ in range(h)] for _ in range(k)], h, p)
        Q = nm([[rng.randrange(p) for _ in range(k)] for _ in range(m)], k, p)
        A1 = tl(P * nm(M1, m, p) * Q)
        A2 = tl(P * nm(M2, m, p) * Q)
        xs = list(range(k + 1))
        ys = [int(nm([[(s * A1[i][j] + A2[i][j]) % p for j in range(k)] for i in range(k)], k, p).det()) for s in xs]
        f = interp(xs, ys, p)
        tops.append(int(nm(A1, k, p).det()))
        g = f if g is None else g.gcd(f)
        if g.degree() == 0 and any(tops):
            return True
    return False


def minrank_lower(M1, M2, p, rng):
    """largest k with a certificate (so min rank >= k)."""
    h = len(M1)
    m = len(M1[0]) if h else 0
    best = 0
    for k in range(1, min(h, m) + 1):
        if minrank_cert(M1, M2, k, p, rng):
            best = k
        else:
            break
    return best


def ub1_block(thS, thT, hs, ht, p, order=None):
    """dim of {phi = (phi_l): phi_{l+1} theta^s_v = theta^t_v phi_l for all v, l} (first-order solution space)."""
    nl = len(hs)
    offs, U = [], 0
    for l in range(nl):
        offs.append(U)
        U += ht[l] * hs[l]
    if U == 0:
        return 0
    nv = len(thS[0]) if nl > 1 else 0
    basis, k = None, U
    for vi in (order if order is not None else range(nv)):
        eqs = []
        for l in range(nl - 1):
            A, B = thS[l][vi], thT[l][vi]
            hs0, hs1, ht0, ht1 = hs[l], hs[l + 1], ht[l], ht[l + 1]
            if hs0 == 0 or ht1 == 0:
                continue
            for i in range(ht1):
                b1 = offs[l + 1] + i * hs1
                Bi = B[i]
                for kk in range(hs0):
                    row = [0] * U
                    nz = False
                    for mm in range(hs1):
                        a = A[mm][kk]
                        if a:
                            row[b1 + mm] = a
                            nz = True
                    for mm in range(ht0):
                        b = Bi[mm]
                        if b:
                            idx = offs[l] + mm * hs0 + kk
                            row[idx] = (row[idx] - b) % p
                            nz = True
                    if nz:
                        eqs.append(row)
        if not eqs:
            continue
        Em = nm(eqs, U, p)
        C = Em if basis is None else Em * basis
        X, nul = C.nullspace()
        if nul == 0:
            return 0
        if nul == k:
            continue
        Xl = tl(X)
        newb = nm([[Xl[i][j] for j in range(nul)] for i in range(k)], nul, p)
        basis = newb if basis is None else basis * newb
        k = nul
    return k


def rel_lattice_Q(S):
    """Q-basis of {n in Q^S : sum_m n_m (m, 1) = 0}."""
    m = len(S)
    A = [[Fraction(S[j][i]) for j in range(m)] for i in range(3)] + [[Fraction(1)] * m]
    # rref over Q
    piv, r = [], 0
    for col in range(m):
        pr = next((i for i in range(r, 4) if A[i][col] != 0), None)
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        iv = 1 / A[r][col]
        A[r] = [x * iv for x in A[r]]
        for i in range(4):
            if i != r and A[i][col] != 0:
                f = A[i][col]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(col)
        r += 1
        if r == 4:
            break
    free = [j for j in range(m) if j not in piv]
    basis = []
    for fj in free:
        v = [Fraction(0)] * m
        v[fj] = Fraction(1)
        for i, pc in enumerate(piv):
            v[pc] = -A[i][fj]
        basis.append(v)
    return basis


def perm_apply(pi, m):
    """alpha: x_i -> lambda_i x_{pi(i)}; exponent m maps to e' with e'[pi(i)] = m[i]."""
    e = [0, 0, 0]
    for i in range(3):
        e[pi[i]] = m[i]
    return tuple(e)


def allowed_perms(S):
    """permutations pi (with a torus factor depending on F) that preserve the very general member up to scalar."""
    Sset = set(S)
    idx = {m: k for k, m in enumerate(S)}
    rel = rel_lattice_Q(S)
    out = []
    for pi in itertools.permutations(range(3)):
        if set(perm_apply(pi, m) for m in S) != Sset:
            continue
        if all(all(v[idx[perm_apply(pi, m)]] == v[idx[m]] for m in S) for v in rel):
            out.append(pi)
    return out


def is_prime(m):
    if m < 2:
        return False
    i = 2
    while i * i <= m:
        if m % i == 0:
            return False
        i += 1
    return True


def lcm(*a):
    r = 1
    for x in a:
        r = r * x // gcd(r, x)
    return r


def pick_primes(L, count=2, start=50000, excl=(31081, 30241, 32003)):
    out = []
    p = start + ((1 - start) % L)
    while len(out) < count:
        if is_prime(p) and p not in excl:
            out.append(p)
        p += L
    return out


def gmax_order(S):
    m0 = S[0]
    vs = [(m[1] - m0[1], m[2] - m0[2]) for m in S[1:]]
    g = 0
    for a in range(len(vs)):
        for b in range(a + 1, len(vs)):
            g = gcd(g, abs(vs[a][0] * vs[b][1] - vs[a][1] * vs[b][0]))
    return g
