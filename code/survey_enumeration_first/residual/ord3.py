"""D1_order3: independent first/second/third-order test for one piece H_i of Y_F = {F(x0,x1,x2) + u^d = 0},
pieces with h^{0,2} = 0 (i in {0,1}), written from the formulas (not from so2.py / f4 code).

H = H_i = A + B at the base point, A = F^2 = span{omega_P = eta(P,0)}, P in S_{d-4-i};
B = span{eta_B(M) = eta(M,1)}, M standard monomials of (R_F)_{2d-4-i}; eta(X,q) = Res(u^i X Omega / G^{q+1}).
Linear family F_t = F + L_t: omega_P(t) = sum_k (-1)^k eta(P L^k, k) =: sum_k W_k(P) (exact Taylor expansion).
Griffiths: eta(sum g_j d_jF, q) = (1/q) eta(sum d_j g_j, q-1); eta(X,1) = eta_B(NF X) + omega_{div h},
X = NF X + sum h_j d_jF.  A flat phi with phi(A_t) in A_t: phi = [[c, beta], [0, phibar]] in A+B,
phi(W_k(P)) = sum_{j<=k} W_{k-j}(c^(j) P)   (c(t) = sum_j c^(j), c^(0) = c, W_0(Q) = omega_Q).
B-component of degree k: constraint; A-component of degree k: determines c^(k).
Variants: V1 = per-arc free c^(j) (j < order), no beta, A-components dropped (= the so2/f4 relaxation at order 2);
          V3 = exact jet condition: beta and A-components imposed (c^(j) per arc determined by them).
Blocks (sigma, Delta) under the diagonal group G = mu_n with weights w (chi(e) = w.e mod n).
Usage: python -I -B ord3.py FAMILY d i p seed narcs order [small]
"""
import json
import os
import random
import sys
import time

import flint

sys.stdout.reconfigure(encoding='utf-8')
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'outputs')


def monos(D):
    if D < 0:
        return []
    return [(a, b, D - a - b) for a in range(D, -1, -1) for b in range(D - a, -1, -1)]


def pmul(P, Q, p):
    out = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = (e1[0] + e2[0], e1[1] + e2[1], e1[2] + e2[2])
            out[e] = (out.get(e, 0) + c1 * c2) % p
    return {e: c for e, c in out.items() if c}


def pder(P, j, p):
    out = {}
    for e, c in P.items():
        if e[j]:
            f = list(e)
            f[j] -= 1
            out[tuple(f)] = (out.get(tuple(f), 0) + c * e[j]) % p
    return {e: c for e, c in out.items() if c}


def family(name, d):
    """basis (list of polys), G = (n, w)."""
    allm = monos(d)
    if name == 'gen':
        return [{e: 1} for e in allm], (1, (0, 0, 0))
    if name in ('pp', 'pm'):   # d-ics even in x1 and x2; pp: also symmetric under x1 <-> x2. G = mu_2 on x2 only
        ev = [e for e in allm if e[1] % 2 == 0 and e[2] % 2 == 0]
        if name == 'pm':
            return [{e: 1} for e in ev], (2, (0, 0, 1))
        seen, basis = set(), []
        for e in ev:
            if e in seen:
                continue
            f = (e[0], e[2], e[1])
            seen |= {e, f}
            basis.append({e: 1} if e == f else {e: 1, f: 1})
        return basis, (2, (0, 0, 1))
    if name.startswith('cyc'):   # cyc<n>_<w0w1w2>[ab|Qs|Qsab]
        n = int(name[3:name.index('_')])
        rest = name[name.index('_') + 1:]
        w = tuple(int(c) for c in rest.split('-')[:3]) if '-' in rest else tuple(int(c) for c in rest[:3])  # Z1: dashed weights
        sub = '' if '-' in rest else rest[3:]
        mons = [e for e in allm if sum(a * b for a, b in zip(w, e)) % n == 0]
        basis = [{e: 1} for e in mons]
        xd1, xd2 = (0, d, 0), (0, 0, d)
        if sub in ('Qs', 'Qsab'):   # Q = s^4 + l s^2 t^2 + t^4 (s = x0^2, t = x1 x2): drop x0^6 x1x2, x0^2 x1^3 x2^3
            basis = [P for P in basis if next(iter(P)) not in ((d - 2, 1, 1), (2, 3, 3))]
        if sub in ('ab', 'Qsab'):
            basis = [P for P in basis if next(iter(P)) not in (xd1, xd2)] + [{xd1: 1, xd2: 1}]
        return basis, (n, w)
    raise ValueError(name)


class Setup:
    def __init__(self, F, d, i0, p):
        self.F, self.d, self.i0, self.p = F, d, i0, p
        self.parts = [pder(F, j, p) for j in range(3)]
        self.D = [d - 4 - i0 + k * d for k in range(4)]   # degrees of eta(., k)
        D0, D1, D2, D3 = self.D
        self.S = [monos(D) for D in self.D]
        self.idx = [{e: k for k, e in enumerate(S)} for S in self.S]
        # --- Red1 on S_{D1}: T = [Jacobian columns | standard monomials]
        H1 = monos(D1 - (d - 1))
        cols = [pmul({m: 1}, self.parts[j], p) for j in range(3) for m in H1]
        nh = len(cols)
        MT = flint.nmod_mat(nh, len(self.S[1]), p)
        for r, P in enumerate(cols):
            for e, c in P.items():
                MT[r, self.idx[1][e]] = c
        R, rk = MT.rref()
        assert rk == nh, "Jacobian map in degree D1 not injective"
        piv = set()
        for r in range(rk):
            piv.add(next(k for k in range(len(self.S[1])) if int(R[r, k])))
        self.std = [k for k in range(len(self.S[1])) if k not in piv]   # indices into S_{D1}
        n1 = len(self.S[1])
        T = flint.nmod_mat(n1, n1, p)
        for cidx, P in enumerate(cols):
            for e, c in P.items():
                T[self.idx[1][e], cidx] = c
        for k, s in enumerate(self.std):
            T[s, nh + k] = 1
        Tinv = T.inv()
        # Div1: S_{D1-d+1}^3 -> S_{D0}
        Div1 = self.divmat(H1, 0)
        self.RA1 = Div1 * self.rows(Tinv, 0, nh)          # |S0| x |S1|
        self.RB1 = self.rows(Tinv, nh, n1)                # |std| x |S1|
        self.nA, self.nB = len(self.S[0]), len(self.std)
        # --- Red2: S_{D2} -> S_{D1}, Red3: S_{D3} -> S_{D2}
        self.Red2 = self.redmat(2, None)
        self.Red3 = self.redmat(3, None)
        self.RA2, self.RB2 = self.RA1 * self.Red2, self.RB1 * self.Red2
        self.RA3, self.RB3 = self.RA2 * self.Red3, self.RB2 * self.Red3

    def rows(self, M, a, b):
        out = flint.nmod_mat(b - a, M.ncols(), self.p)
        for r in range(a, b):
            for c in range(M.ncols()):
                out[r - a, c] = M[r, c]
        return out

    def divmat(self, Hm, ktarget):
        """matrix of (g_0,g_1,g_2) (each over monomials Hm) -> sum_j d_j g_j in S_{D[ktarget]}"""
        p = self.p
        M = flint.nmod_mat(len(self.S[ktarget]), 3 * len(Hm), p)
        for j in range(3):
            for k, m in enumerate(Hm):
                if m[j]:
                    f = list(m)
                    f[j] -= 1
                    M[self.idx[ktarget][tuple(f)], j * len(Hm) + k] = m[j] % p
        return M

    def redmat(self, k, perm_seed):
        """S_{D[k]} -> S_{D[k-1]}: X = sum g_j d_jF (any preimage), X -> div g. perm_seed: shuffle column order
        (different preimage choice) for the consistency check."""
        p, d = self.p, self.d
        Hk = monos(self.D[k] - (d - 1))
        cols = [(j, m) for j in range(3) for m in Hk]
        order = list(range(len(cols)))
        if perm_seed is not None:
            random.Random(perm_seed).shuffle(order)
        nr = len(self.S[k])
        M = flint.nmod_mat(nr, len(cols), p)
        for ci, oc in enumerate(order):
            j, m = cols[oc]
            for e, c in pmul({m: 1}, self.parts[j], p).items():
                M[self.idx[k][e], ci] = c
        R, rk = M.rref()
        assert rk == nr, "R_F nonzero in degree %d (F singular mod p?)" % self.D[k]
        piv = [next(c for c in range(len(cols)) if int(R[r, c])) for r in range(rk)]
        Sub = flint.nmod_mat(nr, nr, p)
        for r in range(nr):
            for a, c in enumerate(piv):
                Sub[r, a] = M[r, c]
        Si = Sub.inv()
        RI = flint.nmod_mat(len(cols), nr, p)    # right inverse, in the ORIGINAL column order
        for a, c in enumerate(piv):
            for r in range(nr):
                RI[order[c], r] = Si[a, r]
        return self.divmat(Hk, k - 1) * RI

    def vec(self, P, k):
        v = flint.nmod_mat(len(self.S[k]), 1, self.p)
        for e, c in P.items():
            v[self.idx[k][e], 0] = c
        return v

    def W(self, P, L, kmax):
        """W_k(P) for k = 1..kmax along L: list of (A-list, B-list)"""
        p = self.p
        out = []
        X = P
        inv2, inv6 = pow(2, p - 2, p), pow(6, p - 2, p)
        for k in range(1, kmax + 1):
            X = pmul(X, L, p)
            v = self.vec(X, k)
            RA, RB, sc = {1: (self.RA1, self.RB1, p - 1), 2: (self.RA2, self.RB2, inv2),
                          3: (self.RA3, self.RB3, (p - inv6) % p)}[k]
            a, b = RA * v, RB * v
            out.append(([int(a[r, 0]) * sc % p for r in range(self.nA)],
                        [int(b[r, 0]) * sc % p for r in range(self.nB)]))
        return out


def nullspace_proj(rows, N, projs, p):
    """rows: list of dicts {col: val}. Returns (kernel dim, [rank of projection to coords [0, m)) for m in projs])"""
    if not rows:
        rows = [{}]
    M = flint.nmod_mat(len(rows), N, p)
    for r, row in enumerate(rows):
        for c, v in row.items():
            if v % p:
                M[r, c] = v % p
    X, nul = M.nullspace()
    out = []
    for m in projs:
        if nul == 0 or m == 0:
            out.append(0)
            continue
        Y = flint.nmod_mat(nul, m, p)
        for a in range(nul):
            for r in range(m):
                Y[a, r] = X[r, a]
        out.append(Y.rank())
    return nul, out, (X, nul)


def cbasis(X, nul, nc, p):
    Y = flint.nmod_mat(nul, nc, p)
    for a in range(nul):
        for r in range(nc):
            Y[a, r] = X[r, a]
    R, rk = Y.rref()
    return [[int(R[a, r]) for r in range(nc)] for a in range(rk)]


def main():
    name, d, i0, p, seed, narcs, order = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), \
        int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
    small = len(sys.argv) > 8 and sys.argv[8] == 'small'
    t0 = time.time()
    rng = random.Random(777000 + 13 * seed + p)
    basis, (n, w) = family(name, d)
    F = {}
    for P in basis:
        c = rng.randrange(1, 40) if small else rng.randrange(1, p)
        for e, x in P.items():
            F[e] = (F.get(e, 0) + c * x) % p
    S = Setup(F, d, i0, p)
    # consistency check: Griffiths reduction independent of the preimage choice
    for k in (2, 3):
        alt = S.redmat(k, 99 + k)
        Xs = [{e: rng.randrange(p) for e in S.S[k]} for _ in range(3)]
        for X in Xs:
            v = S.vec(X, k)
            if k == 2:
                r1, r2 = S.Red2 * v, alt * v
                a1, b1, a2, b2 = S.RA1 * r1, S.RB1 * r1, S.RA1 * r2, S.RB1 * r2
            else:
                r1, r2 = S.Red2 * (S.Red3 * v), S.Red2 * (alt * v)
                a1, b1, a2, b2 = S.RA1 * r1, S.RB1 * r1, S.RA1 * r2, S.RB1 * r2
            assert a1 == a2 and b1 == b2, "reduction depends on preimage choice (k=%d)" % k
    print("%s d=%d i=%d p=%d seed=%d narcs=%d order=%d small=%s | dimA=%d dimB=%d |basis|=%d G=mu_%d w=%s "
          "(reduction consistency check passed; setup %.1fs)" % (name, d, i0, p, seed, narcs, order, small, S.nA,
                                                                 S.nB, len(basis), n, w, time.time() - t0))
    chi = lambda e: sum(a * b for a, b in zip(w, e)) % n
    chA = [chi(e) for e in S.S[0]]
    chB = [chi(S.S[1][s]) for s in S.std]

    def rand_dir():
        L = {}
        for P in basis:
            c = rng.randrange(p)
            for e, x in P.items():
                L[e] = (L.get(e, 0) + c * x) % p
        return {e: c for e, c in L.items() if c}
    if len(basis) > 20:
        tangents = [rand_dir() for _ in range(12)]
    else:
        tangents = [dict(P) for P in basis]
    arcs = [rand_dir() for _ in range(narcs)]
    Amons = [{e: 1} for e in S.S[0]]
    WT = [[S.W(P, T, 1)[0] for P in Amons] for T in tangents]          # WT[t][P] = (A, B) order 1
    WA = [[S.W(P, L, order) for P in Amons] for L in arcs]              # WA[k][P][ord-1] = (A, B)
    if 'sab3' in sys.argv[8:] and order >= 3:   # deliberate sabotage (control): W_3 -> 2 W_3
        WA = [[[w if o < 2 else ([2 * x % p for x in w[0]], [2 * x % p for x in w[1]]) for o, w in enumerate(WP)]
               for WP in Wk] for Wk in WA]
        print("SABOTAGE: W_3 doubled")
    print("W computed (%.1fs)" % (time.time() - t0))
    res = {"family": name, "d": d, "i": i0, "p": p, "seed": seed, "narcs": narcs, "order": order, "small": small,
           "dimA": S.nA, "dimB": S.nB, "G": [n, list(w)], "F": {str(e): c for e, c in F.items()}, "blocks": []}

    def system(sA, tA, sB, tB, variant, ordmax):
        nc, nph = len(tA) * len(sA), len(tB) * len(sB)
        nbe = len(tA) * len(sB) if (variant == 'V3' and ordmax >= 2) else 0
        nj = max(ordmax - 1, 0)
        N = nc + nph + nbe + nc * nj * narcs
        ic = lambda t, s: t * len(sA) + s
        iph = lambda t, s: nc + t * len(sB) + s
        ibe = lambda t, s: nc + nph + t * len(sB) + s
        icd = lambda k, j, t, s: nc + nph + nbe + ((k * nj + (j - 1)) * nc) + t * len(sA) + s

        def cvar(k, j, t, s):
            return ic(t, s) if j == 0 else icd(k, j, t, s)
        rows = []
        for Tw in WT:   # order 1 along the basis tangents
            for s, P in enumerate(sA):
                for tbi, tb in enumerate(tB):
                    row = {}
                    for sbi, sb in enumerate(sB):
                        if Tw[P][1][sb]:
                            row[iph(tbi, sbi)] = Tw[P][1][sb]
                    for t, Q in enumerate(tA):
                        x = Tw[Q][1][tb]
                        if x:
                            row[ic(t, s)] = (row.get(ic(t, s), 0) - x) % p
                    if row:
                        rows.append(row)
        for k in range(narcs):
            Wk = WA[k]
            for o in range(1, ordmax + 1):
                for s, P in enumerate(sA):
                    WP = Wk[P][o - 1]
                    for tbi, tb in enumerate(tB):     # B-component of degree o
                        row = {}
                        for sbi, sb in enumerate(sB):
                            if WP[1][sb]:
                                row[iph(tbi, sbi)] = WP[1][sb]
                        for j in range(o):
                            for t, Q in enumerate(tA):
                                x = Wk[Q][o - j - 1][1][tb]
                                if x:
                                    vv = cvar(k, j, t, s)
                                    row[vv] = (row.get(vv, 0) - x) % p
                        if row:
                            rows.append(row)
                    if variant == 'V3' and o <= ordmax - 1:   # A-component of degree o (defines c^(o))
                        for tai, ta in enumerate(tA):
                            row = {icd(k, o, tai, s): p - 1}
                            for s2, P2 in enumerate(sA):
                                x = WP[0][P2]
                                if x:
                                    row[ic(tai, s2)] = (row.get(ic(tai, s2), 0) + x) % p
                            for sbi, sb in enumerate(sB):
                                if WP[1][sb]:
                                    row[ibe(tai, sbi)] = (row.get(ibe(tai, sbi), 0) + WP[1][sb]) % p
                            for j in range(o):
                                for t, Q in enumerate(tA):
                                    x = Wk[Q][o - j - 1][0][ta]
                                    if x:
                                        vv = cvar(k, j, t, s)
                                        row[vv] = (row.get(vv, 0) - x) % p
                            rows.append(row)
        nul, (pc, pcp), (X, nl) = nullspace_proj(rows, N, (nc, nc + nph), p)
        return nul, pc, pcp, (X, nl, nc)

    variants = [('o1', 'V1', 1)] + [('V1o%d' % o, 'V1', o) for o in range(2, order + 1)] + \
               [('V3o%d' % o, 'V3', o) for o in range(2, order + 1)]
    for sg in range(n):
        for De in range(n):
            sA = [k for k in range(S.nA) if chA[k] == sg]
            tA = [k for k in range(S.nA) if chA[k] == (sg + De) % n]
            sB = [k for k in range(S.nB) if chB[k] == sg]
            tB = [k for k in range(S.nB) if chB[k] == (sg + De) % n]
            if not sA or not tA:
                continue
            blk = {"sigma": sg, "Delta": De, "c_shape": [len(tA), len(sA)], "phibar_shape": [len(tB), len(sB)]}
            last = None
            for lab, var, o in variants:
                nul, pc, pcp, last_ = system(sA, tA, sB, tB, var, o)
                blk[lab] = {"c": pc, "c_phibar": pcp, "full": nul}
                if lab == 'V3o%d' % order or (order == 1 and lab == 'o1'):
                    last = last_
            X, nl, nc = last
            if nc <= 6 and nl:
                blk["c_basis_last"] = cbasis(X, nl, nc, p)
            blk["tA_mons"] = [list(S.S[0][k]) for k in tA]
            blk["sA_mons"] = [list(S.S[0][k]) for k in sA]
            res["blocks"].append(blk)
            if any(blk[lab]["c"] for lab, _, _ in variants):
                print("sigma=%d Delta=%d c %dx%d phibar %dx%d | c-proj: %s | (c,phibar)-proj: %s | full: %s%s" % (
                    sg, De, len(tA), len(sA), len(tB), len(sB),
                    " ".join("%s=%d" % (lab, blk[lab]["c"]) for lab, _, _ in variants),
                    " ".join("%s=%d" % (lab, blk[lab]["c_phibar"]) for lab, _, _ in variants),
                    " ".join("%s=%d" % (lab, blk[lab]["full"]) for lab, _, _ in variants),
                    ("  c-basis(last)=%s" % blk.get("c_basis_last")) if blk.get("c_basis_last") else ""))
    print("done %.1fs" % (time.time() - t0))
    fn = os.path.join(OUT, "ord3_%s_d%d_i%d_p%d_s%d_a%d_o%d%s.json" % (name, d, i0, p, seed, narcs, order,
                                                                       "_small" if small else ""))
    with open(fn, "w", newline="\n") as fh:
        fh.write(json.dumps(res, indent=1) + "\n")


if __name__ == '__main__':
    main()
