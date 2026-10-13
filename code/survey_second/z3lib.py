# z3lib.py -- agent Z3, independent toolkit (written from TOOLS.txt / S2 / THEOREM_U lemmas only).
# Jacobian ring mod p by (degree, G-character), Griffiths reduction, pieces E_(a,chi), Higgs/multiplication maps,
# first-order (Hom-)commutants, flat Taylor coefficients, J-system, hull ranks, L4' pencil test,
# Galois orbits, normaliser permutations, curve planes, eigen-split at a symmetric point.
import sys, math, random, itertools
import numpy as np
from flint import nmod_mat, nmod_poly

INF = 10**9

# ------------------------------------------------------------------ basics
def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    dd, s = n - 1, 0
    while dd % 2 == 0:
        dd //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, dd, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def lcm(*xs):
    r = 1
    for x in xs:
        r = r * x // math.gcd(r, x)
    return r

def pick_primes(d, n, count=2, lo=100003, avoid=(90001, 90217)):
    M = lcm(8, d, n)
    p = lo + ((1 - lo) % M)
    out = []
    while len(out) < count:
        if is_prime(p) and p not in avoid and not (49000 <= p <= 51000):
            out.append(p)
        p += M
    return out

_MON = {}
def monos(deg):
    if deg < 0:
        return []
    if deg not in _MON:
        _MON[deg] = [(i, j, deg - i - j) for i in range(deg, -1, -1) for j in range(deg - i, -1, -1)]
    return _MON[deg]

def parse_family(name):
    s = name[3:]
    ns, ws = s.split('_')
    n = int(ns)
    w = tuple(int(x) for x in ws.split('-')) if '-' in ws else tuple(int(ch) for ch in ws)
    return n, w

def wch(b, w, n):
    return (w[0] * b[0] + w[1] * b[1] + w[2] * b[2]) % n

def family_S(d, n, w):
    return [m for m in monos(d) if wch(m, w, n) == 0]

def addm(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])

def perm_exp(m, pi):
    # (pi m)_{pi(i)} = m_i
    r = [0, 0, 0]
    for i in range(3):
        r[pi[i]] = m[i]
    return tuple(r)

# ------------------------------------------------------------------ flint <-> numpy
def F2N(M):
    r, c = M.nrows(), M.ncols()
    if r == 0 or c == 0:
        return np.zeros((r, c), dtype=np.int64)
    try:
        ent = M.entries()
    except AttributeError:
        ent = [M[i, j] for i in range(r) for j in range(c)]
    return np.array([int(x) for x in ent], dtype=np.int64).reshape(r, c)

def N2F(A, p):
    r, c = A.shape
    return nmod_mat(r, c, [int(x) for x in (A % p).ravel().tolist()], p)

def rank_mod(A, p):
    if A.shape[0] == 0 or A.shape[1] == 0:
        return 0
    return N2F(A, p).rank()

def nullity_mod(A, p):
    return A.shape[1] - rank_mod(A, p)

def nullspace_mod(A, p):
    r, c = A.shape
    if c == 0:
        return np.zeros((0, 0), dtype=np.int64)
    if r == 0:
        return np.eye(c, dtype=np.int64)
    X, nul = N2F(A, p).nullspace()
    return F2N(X)[:, :nul]

def inv_mod(A, p):
    if A.shape[0] == 0:
        return A.copy()
    return F2N(N2F(A, p).inv())

def det_mod(A, p):
    if A.shape[0] == 0:
        return 1
    return int(N2F(A, p).det())

def rref_pivots(A, p):
    if A.shape[0] == 0 or A.shape[1] == 0:
        return [], 0
    R, rk = N2F(A, p).rref()
    Rn = F2N(R)[:rk]
    piv = [int(np.argmax(Rn[i] != 0)) for i in range(rk)]
    return piv, rk

def poly_mul(A, B, p):
    out = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            k = (ma[0] + mb[0], ma[1] + mb[1], ma[2] + mb[2])
            out[k] = (out.get(k, 0) + ca * cb) % p
    return {k: v for k, v in out.items() if v}

# ------------------------------------------------------------------ Jacobian ring
class Ring:
    """R_F = k[x0,x1,x2]/(dF) over F_p, graded by (degree e, G-character c = w.m mod n).
    red(e,c): standard monomials (basis of R_{e,c}), NF (normal form of each monomial), Low (Griffiths lowering:
    x^mu = NF(mu) + sum_j g_j d_jF  ->  sum_j d_j g_j, degree e-d)."""
    def __init__(self, d, n, w, F, p):
        self.d, self.n, self.w, self.p = d, n, tuple(w), p
        self.F = {m: c % p for m, c in F.items() if c % p}
        self.S = family_S(d, n, w)
        self.dF = []
        for j in range(3):
            D = {}
            for m, c in self.F.items():
                if m[j] > 0:
                    mm = list(m); mm[j] -= 1
                    D[tuple(mm)] = (c * m[j]) % p
            self.dF.append(D)
        self._mons = {}
        self._red = {}

    def mons(self, e, c):
        key = (e, c)
        r = self._mons.get(key)
        if r is None:
            L = [m for m in monos(e) if wch(m, self.w, self.n) == c]
            r = (L, {m: i for i, m in enumerate(L)})
            self._mons[key] = r
        return r

    def red(self, e, c):
        key = (e, c)
        if key in self._red:
            return self._red[key]
        d, n, w, p = self.d, self.n, self.w, self.p
        ML, MI = self.mons(e, c)
        N = len(ML)
        LL, LI = self.mons(e - d, c)
        gens = []
        if e >= d - 1:
            for j in range(3):
                for nu in self.mons(e - d + 1, (c + w[j]) % n)[0]:
                    gens.append((j, nu))
        if N == 0:
            dat = dict(std=[], sidx={}, NF=np.zeros((0, 0), np.int64), Low=np.zeros((len(LL), 0), np.int64), N=0, MI=MI)
        elif not gens:
            dat = dict(std=list(ML), sidx={m: i for i, m in enumerate(ML)}, NF=np.eye(N, dtype=np.int64),
                       Low=np.zeros((len(LL), N), np.int64), N=N, MI=MI)
        else:
            K = len(gens)
            M = np.zeros((N, K), np.int64)
            for col, (j, nu) in enumerate(gens):
                for mu, cf in self.dF[j].items():
                    r = MI[(nu[0] + mu[0], nu[1] + mu[1], nu[2] + mu[2])]
                    M[r, col] = (M[r, col] + cf) % p
            pivm, rk = rref_pivots(M.T.copy(), p)
            pset = set(pivm)
            std_idx = [i for i in range(N) if i not in pset]
            pivg, rk2 = rref_pivots(M, p)
            assert rk2 == rk
            A = np.zeros((N, N), np.int64)
            A[:, :rk] = M[:, pivg]
            for t, i in enumerate(std_idx):
                A[i, rk + t] = 1
            Ainv = inv_mod(A, p)
            NF = Ainv[rk:, :].copy()
            G = Ainv[:rk, :]
            Lg = np.zeros((len(LL), rk), np.int64)
            for t, gi in enumerate(pivg):
                j, nu = gens[gi]
                if nu[j] > 0:
                    mm = list(nu); mm[j] -= 1
                    Lg[LI[tuple(mm)], t] = nu[j] % p
            Low = (Lg @ G) % p
            std = [ML[i] for i in std_idx]
            dat = dict(std=std, sidx={m: i for i, m in enumerate(std)}, NF=NF, Low=Low, N=N, MI=MI)
        self._red[key] = dat
        return dat

    def hilb(self, e, c):
        return len(self.red(e, c)['std'])

    def koszul(self, e, c):
        tot = 0
        for k in range(4):
            for J in itertools.combinations(range(3), k):
                e2 = e - (self.d - 1) * k
                c2 = (c + sum(self.w[j] for j in J)) % self.n
                tot += (-1) ** k * len(self.mons(e2, c2)[0])
        return tot

    def koszul_check(self):
        bad = []
        for e in range(0, 3 * self.d - 4):
            for c in range(self.n):
                if self.hilb(e, c) != self.koszul(e, c):
                    bad.append((e, c))
        return bad

# ------------------------------------------------------------------ pieces
class Piece:
    """E_(a,chi): residues Res(x^b u^{a-1} Omega / G^{q+1}), |b| = D0 + q d, chi = w.b + sum(w) mod n."""
    def __init__(self, R, a, chi):
        self.R, self.a, self.chi = R, a, chi
        d, n, w = R.d, R.n, R.w
        self.cc = (chi - sum(w)) % n
        self.D0 = d - 3 - a
        self.deg = [self.D0 + q * d for q in range(3)]
        self.lev = [R.red(e, self.cc)['std'] if e >= 0 else [] for e in self.deg]
        self.h = tuple(len(L) for L in self.lev)
        self.dim = sum(self.h)
        self.off = (0, self.h[0], self.h[0] + self.h[1])
        self._mu = None

    def mu(self):
        """list over m in S of [mu_m: level0->level1, level1->level2] (normal-form matrices)."""
        if self._mu is None:
            R = self.R
            dats = [R.red(e, self.cc) if e >= 0 else None for e in self.deg]
            out = []
            for m in R.S:
                pair = []
                for q in (0, 1):
                    if self.h[q] and self.h[q + 1]:
                        tgt = dats[q + 1]
                        cols = [tgt['MI'][addm(b, m)] for b in self.lev[q]]
                        pair.append(tgt['NF'][:, cols].copy())
                    else:
                        pair.append(np.zeros((self.h[q + 1], self.h[q]), np.int64))
                out.append(pair)
            self._mu = out
        return self._mu

# ------------------------------------------------------------------ commutant systems
def hom_system(hx, hy, mus_x, mus_y, p):
    """Unknowns X_q: E_x^q -> E_y^q (column-major vec); equations X_{q+1} mu^x_q = mu^y_q X_q."""
    nunk = [hy[q] * hx[q] for q in range(3)]
    offs = [0, nunk[0], nunk[0] + nunk[1]]
    N = sum(nunk)
    blocks = []
    for mx, my in zip(mus_x, mus_y):
        for q in (0, 1):
            rows = hy[q + 1] * hx[q]
            if rows == 0:
                continue
            M = np.zeros((rows, N), np.int64)
            if nunk[q + 1]:
                M[:, offs[q + 1]:offs[q + 1] + nunk[q + 1]] = np.kron(mx[q].T, np.eye(hy[q + 1], dtype=np.int64))
            if nunk[q]:
                M[:, offs[q]:offs[q] + nunk[q]] = (-np.kron(np.eye(hx[q], dtype=np.int64), my[q])) % p
            blocks.append(M)
    if not blocks:
        return np.zeros((0, N), np.int64)
    return np.vstack(blocks) % p

def hom_commutant_dim(Px, Py, p):
    A = hom_system(Px.h, Py.h, Px.mu(), Py.mu(), p)
    return A.shape[1] - rank_mod(A, p)

# ------------------------------------------------------------------ Griffiths reduction and Taylor coefficients
def reduce_mat(R, V, e, c, q):
    """class of Res(V u^{a-1} Omega / G^q) (columns of V over mons(e,c)) -> {level: coords over std}."""
    p, d = R.p, R.d
    out = {}
    while True:
        dat = R.red(e, c)
        nf = (dat['NF'] @ V) % p if dat['N'] else np.zeros((0, V.shape[1]), np.int64)
        if q - 1 <= 2:
            out[q - 1] = nf
        elif nf.size and nf.any():
            raise RuntimeError('nonzero normal form above the socle')
        if q == 1:
            break
        V = ((dat['Low'] @ V) % p) * pow(q - 1, p - 2, p) % p
        e -= d
        q -= 1
    return out

def place(P, red, r):
    W = np.zeros((P.dim, r), np.int64)
    for lev, M in red.items():
        if 0 <= lev <= 2 and P.h[lev]:
            W[P.off[lev]:P.off[lev] + P.h[lev], :] = M
    return W

class Arc:
    def __init__(self, L, K, p):
        self.L = L
        self.pows = [{(0, 0, 0): 1}]
        for _ in range(K):
            self.pows.append(poly_mul(self.pows[-1], L, p))

def taylor_W(R, P, arc, K):
    """W_k(P) = (-1)^k [Res(P L^k Omega/G^{k+1})], k = 0..K, for the top basis P of piece P (columns)."""
    p, d = R.p, R.d
    r = P.h[0]
    out = []
    for k in range(K + 1):
        e = P.D0 + k * d
        ML, MI = R.mons(e, P.cc)
        V = np.zeros((len(ML), r), np.int64)
        for col, b in enumerate(P.lev[0]):
            for mu, cf in arc.pows[k].items():
                i = MI[addm(b, mu)]
                V[i, col] = (V[i, col] + cf) % p
        W = place(P, reduce_mat(R, V, e, P.cc, k + 1), r)
        if k % 2:
            W = (-W) % p
        out.append(W)
    return out

def taylor_general(R, P, deltas, K):
    """Arc F_t = F + sum_{j>=1} t^j deltas[j-1] (deltas G-invariant of degree d). W_k = sum_r (-1)^r [t^k] delta^r / G^{r+1}."""
    p, d = R.p, R.d
    Dp = [[{(0, 0, 0): 1}] + [{} for _ in range(K)]]
    for rr in range(1, K + 1):
        prev = Dp[-1]
        cur = [{} for _ in range(K + 1)]
        for k1 in range(K + 1):
            if not prev[k1]:
                continue
            for j in range(1, K + 1 - k1):
                if j - 1 >= len(deltas):
                    break
                prod = poly_mul(prev[k1], deltas[j - 1], p)
                for m, c in prod.items():
                    cur[k1 + j][m] = (cur[k1 + j].get(m, 0) + c) % p
        Dp.append(cur)
    r = P.h[0]
    out = []
    for k in range(K + 1):
        W = np.zeros((P.dim, r), np.int64)
        for rr in range(0, k + 1):
            poly = Dp[rr][k]
            if not poly:
                continue
            e = P.D0 + rr * d
            ML, MI = R.mons(e, P.cc)
            V = np.zeros((len(ML), r), np.int64)
            for col, b in enumerate(P.lev[0]):
                for mu, cf in poly.items():
                    i = MI[addm(b, mu)]
                    V[i, col] = (V[i, col] + cf) % p
            Wr = place(P, reduce_mat(R, V, e, P.cc, rr + 1), r)
            W = (W + (Wr if rr % 2 == 0 else -Wr)) % p
        out.append(W)
    return out

def hull_rank(Wlist, p):
    cols = [W for Wa in Wlist for W in Wa if W.shape[1]]
    if not cols:
        return 0
    return rank_mod(np.hstack(cols), p)

def jtest(Px, Py, Wx, Wy, K, p):
    """J_K kernel dimension (T2 Def. 2.2). Wx[a][k]: dim_x x r_x ; Wy[a][k]: dim_y x r_y."""
    dx, dy, rx, ry = Px.dim, Py.dim, Px.h[0], Py.h[0]
    A = len(Wx)
    nphi, nc = dx * dy, rx * ry
    ncols = nphi + nc * (1 + A * K)
    rows = []
    Iy = np.eye(dy, dtype=np.int64)
    Irx = np.eye(rx, dtype=np.int64)
    for a in range(A):
        for k in range(K + 1):
            M = np.zeros((dy * rx, ncols), np.int64)
            M[:, :nphi] = np.kron(Wx[a][k].T, Iy)
            if nc:
                for m in range(0, k + 1):
                    blk = (-np.kron(Irx, Wy[a][k - m])) % p
                    col0 = nphi if m == 0 else nphi + nc * (1 + a * K + (m - 1))
                    M[:, col0:col0 + nc] = (M[:, col0:col0 + nc] + blk) % p
            rows.append(M)
    Mall = np.vstack(rows) % p
    return ncols - rank_mod(Mall, p)

# ------------------------------------------------------------------ L4' pencil test (h20 = 2)
def interp_mod(xs, ys, p):
    n = len(xs)
    V = np.array([[pow(x, j, p) for j in range(n)] for x in xs], dtype=np.int64)
    Vi = inv_mod(V, p)
    return [int(v) for v in (Vi @ np.array(ys, dtype=np.int64).reshape(n, 1) % p).ravel()]

def l4p_test(P, rng, trials=4):
    """True if min over omega in P^1(F_p-bar) of rank(v -> mu_v omega) > h11/2 (emptiness via gcd of random
    combinations of (k+1)-minors, incl. the point at infinity)."""
    R = P.R
    p = R.p
    if P.h[0] != 2:
        return None
    h11 = P.h[1]
    k1 = h11 // 2 + 1
    S = R.S
    if k1 > min(h11, len(S)):
        return False
    mus = P.mu()
    M1 = np.column_stack([mus[i][0][:, 0] for i in range(len(S))]) % p
    M2 = np.column_stack([mus[i][0][:, 1] for i in range(len(S))]) % p
    g = None
    lead_nonzero = False
    for _ in range(trials):
        A = np.array([[rng.randrange(p) for _ in range(h11)] for _ in range(k1)], dtype=np.int64)
        B = np.array([[rng.randrange(p) for _ in range(k1)] for _ in range(len(S))], dtype=np.int64)
        X1 = ((A @ M1) % p @ B) % p
        X2 = ((A @ M2) % p @ B) % p
        xs = list(range(k1 + 1))
        ys = [det_mod((s * X1 + X2) % p, p) for s in xs]
        f = nmod_poly(interp_mod(xs, ys, p), p)
        if det_mod(X1, p) % p:
            lead_nonzero = True
        g = f if g is None else g.gcd(f)
    return (g is not None) and g.degree() == 0 and lead_nonzero

# ------------------------------------------------------------------ group-theoretic data
def units(L):
    return [t for t in range(1, max(L, 2)) if math.gcd(t, L) == 1] if L > 1 else [1]

def gal(t, x, d, n):
    return ((t * x[0]) % d, (t * x[1]) % n)

def gmax_check(d, n, w, S):
    """|H| = [Z^3 : span S] (diagonal automorphisms of all Y_F, mod scalars) vs |image of Gamma|."""
    M = np.array(S, dtype=object)
    g = 0
    for rows in itertools.combinations(range(len(S)), 3):
        a = [list(S[r]) for r in rows]
        det = (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
               + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))
        g = math.gcd(g, abs(det))
    Lc = lcm(d, n)
    elems = set()
    for j in range(d):
        for k in range(n):
            v = tuple(((-j * (Lc // d)) + k * w[i] * (Lc // n)) % Lc for i in range(3))
            elems.add(v)
    return g, len(elems)

def u_swap_flags(d, S):
    out = []
    for i in range(3):
        if all(m[i] in (0, d) for m in S):
            out.append(i)
    return out

def curve_planes(d, n, w, S):
    planes = []
    for i in range(3):
        Si = [m for m in S if m[i] == 0]
        if len(Si) != 1:
            continue
        M = Si[0]
        j, k = [t for t in range(3) if t != i]
        al, be = M[j], M[k]
        g = math.gcd(d, math.gcd(al, be))
        if g < 2:
            continue
        e = d // g
        wMg = (w[j] * al + w[k] * be) // g
        planes.append(dict(i=i, M=M, g=g, e=e, wMg=wMg))
    return planes

def curve_delta(planes, a, chi, d, n):
    dl = 0
    for pl in planes:
        e = pl['e']
        if a % e == 0 and a % d != 0:
            if (-(a // e) * pl['wMg']) % n == chi % n:
                dl += 1
    return dl

def relation_basis(S):
    import sympy
    M = sympy.Matrix([[m[i] for m in S] for i in range(3)])
    return [list(v) for v in M.nullspace()]

def valid_perms(d, n, w, S):
    """pi != id with (i) pi(S) = S, (ii) pi trivial on the relation lattice (over Q), (iii) well-defined on pieces."""
    Sset = set(S)
    idx = {m: i for i, m in enumerate(S)}
    rel = relation_basis(S)
    out = []
    for pi in itertools.permutations(range(3)):
        if pi == (0, 1, 2):
            continue
        if any(perm_exp(m, pi) not in Sset for m in S):
            continue
        ok2 = True
        for v in rel:
            for m in S:
                if v[idx[perm_exp(m, pi)]] != v[idx[m]]:
                    ok2 = False
        cmap = {}
        ok3 = True
        sw = sum(w)
        for a in range(1, d):
            D0 = d - 3 - a
            for q in range(3):
                e = D0 + q * d
                for b in monos(e):
                    chi = (wch(b, w, n) + sw) % n
                    chi2 = (wch(perm_exp(b, pi), w, n) + sw) % n
                    key = (a, chi)
                    if cmap.setdefault(key, chi2) != chi2:
                        ok3 = False
        order = 2 if all(pi[pi[i]] == i for i in range(3)) else 3
        out.append(dict(pi=pi, i=True, ii=ok2, iii=ok3, order=order,
                        cmap={k: v for k, v in cmap.items()} if ok3 else None))
    return out

# ------------------------------------------------------------------ eigen-split at a symmetric point
def split_piece(Rs, P, pi, p):
    """Pi = substitution by pi on Gr E_x (pi x = x). Returns dict with eigenbases, restricted mus, checks."""
    Pis = []
    for q in range(3):
        if P.h[q] == 0:
            Pis.append(np.zeros((0, 0), np.int64)); continue
        dat = Rs.red(P.deg[q], P.cc)
        cols = [dat['MI'][perm_exp(b, pi)] for b in P.lev[q]]
        Pis.append(dat['NF'][:, cols].copy())
    sq_ok = all(((Pi @ Pi) % p == np.eye(Pi.shape[0], dtype=np.int64)).all() for Pi in Pis)
    comm_ok = True
    for mu in P.mu():
        for q in (0, 1):
            if P.h[q] and P.h[q + 1]:
                if (((Pis[q + 1] @ mu[q]) - (mu[q] @ Pis[q])) % p).any():
                    comm_ok = False
    bases = {}
    for sg in (1, -1):
        bases[sg] = [nullspace_mod((Pis[q] - sg * np.eye(P.h[q], dtype=np.int64)) % p, p) if P.h[q]
                     else np.zeros((0, 0), np.int64) for q in range(3)]
    Tinv = []
    for q in range(3):
        if P.h[q]:
            T = np.hstack([bases[1][q], bases[-1][q]])
            Tinv.append(inv_mod(T, p))
        else:
            Tinv.append(np.zeros((0, 0), np.int64))
    parts = {}
    for sg in (1, -1):
        k = tuple(bases[sg][q].shape[1] if P.h[q] else 0 for q in range(3))
        start = [0 if sg == 1 else (bases[1][q].shape[1] if P.h[q] else 0) for q in range(3)]
        mus = []
        for mu in P.mu():
            pair = []
            for q in (0, 1):
                if k[q] and k[q + 1]:
                    X = (Tinv[q + 1] @ ((mu[q] @ bases[sg][q]) % p)) % p
                    pair.append(X[start[q + 1]:start[q + 1] + k[q + 1], :].copy())
                else:
                    pair.append(np.zeros((k[q + 1], k[q]), np.int64))
            mus.append(pair)
        parts[sg] = dict(h=k, mus=mus)
    return dict(sq_ok=sq_ok, comm_ok=comm_ok, parts=parts)
