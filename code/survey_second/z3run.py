# z3run.py -- agent Z3: per-family Y-part closure (lb/ub per block orbit) from own toolkit z3lib.
# Usage: python -I -B z3run.py [nproc] [family:d ...]
import sys, os, json, time, random, math, traceback
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from z3lib import *

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'outputs')
K = 4
NARCS = 6
MAX_HOM = 4
MAX_J = 3
MAX_SPLIT = 4

FAMILIES = [('cyc4_013', 6), ('cyc5_012', 6), ('cyc5_014', 6), ('cyc6_013', 6), ('cyc8_013', 6),
            ('cyc7_013', 7), ('cyc9_013', 7),
            ('cyc3_012', 8), ('cyc4_012', 8), ('cyc6_015', 8), ('cyc7_012', 8), ('cyc7_013', 8), ('cyc7_016', 8),
            ('cyc8_013', 8), ('cyc8_014', 8), ('cyc12_015', 8), ('cyc16_019', 8),
            ('cyc7_016', 9), ('cyc8_016', 9), ('cyc12_014', 9),
            ('cyc4_013', 10), ('cyc5_012', 10), ('cyc6_013', 10), ('cyc8_017', 10), ('cyc9_012', 10), ('cyc9_018', 10),
            ('cyc10_012', 10), ('cyc10_013', 10), ('cyc15_016', 10), ('cyc16_017', 10), ('cyc20_0-1-11', 10)]
CONTROLS = [('cyc12_017', 6), ('cyc14_018', 7)]


class PrimeData:
    def __init__(self, name, d, p):
        n, w = parse_family(name)
        S = family_S(d, n, w)
        self.p = p
        for tries in range(40):
            rng = random.Random(f'Z3run|{name}|{d}|{p}|{tries}')
            F = {m: rng.randrange(1, p) for m in S}
            R = Ring(d, n, w, F, p)
            if not R.koszul_check():
                break
        else:
            raise RuntimeError('no smooth point')
        self.R, self.tries, self.rng = R, tries, rng
        self.pieces = {}
        for a in range(1, d):
            for chi in range(n):
                x = Piece(R, a, chi)
                if x.dim:
                    self.pieces[(a, chi)] = x
        self.arcs = [Arc({m: rng.randrange(p) for m in S}, K, p) for _ in range(NARCS)]
        self._c, self._h, self._W, self._hull, self._J = {}, {}, {}, {}, {}

    def comm(self, x):
        if x not in self._c:
            self._c[x] = hom_commutant_dim(self.pieces[x], self.pieces[x], self.p)
        return self._c[x]

    def hom(self, x, y):
        if (x, y) not in self._h:
            self._h[(x, y)] = hom_commutant_dim(self.pieces[x], self.pieces[y], self.p)
        return self._h[(x, y)]

    def W(self, x):
        if x not in self._W:
            self._W[x] = [taylor_W(self.R, self.pieces[x], A, K) for A in self.arcs]
        return self._W[x]

    def hull(self, x):
        if x not in self._hull:
            self._hull[x] = hull_rank(self.W(x), self.p) if self.pieces[x].h[0] else 0
        return self._hull[x]

    def J(self, x, y):
        if (x, y) not in self._J:
            self._J[(x, y)] = jtest(self.pieces[x], self.pieces[y], self.W(x), self.W(y), K, self.p)
        return self._J[(x, y)]


class SplitData:
    """Eigen-split under an involution pi at a random pi-symmetric point (T3)."""
    def __init__(self, name, d, p, v, keys):
        n, w = parse_family(name)
        S = family_S(d, n, w)
        pi = v['pi']
        self.p, self.pi = p, pi
        for tries in range(40):
            rng = random.Random(f'Z3split|{name}|{d}|{p}|{pi}|{tries}')
            F = {}
            for m in S:
                if m in F:
                    continue
                c = rng.randrange(1, p)
                F[m] = c
                F[perm_exp(m, pi)] = c
            R = Ring(d, n, w, F, p)
            if not R.koszul_check():
                break
        else:
            raise RuntimeError('no smooth symmetric point')
        self.R, self.tries = R, tries
        self.stable = {x for x in keys if v['cmap'][x] == x[1]}
        self.pieces = {x: Piece(R, x[0], x[1]) for x in keys}
        self.sp = {}
        self.checks = []
        for x in sorted(self.stable):
            s = split_piece(R, self.pieces[x], pi, p)
            self.sp[x] = s
            self.checks.append((x, s['sq_ok'], s['comm_ok']))
        self._sub = {}

    def parts(self, x):
        return [1, -1] if x in self.stable else ['all']

    def ph(self, x, a):
        return self.sp[x]['parts'][a]['h'] if a != 'all' else self.pieces[x].h

    def pmus(self, x, a):
        return self.sp[x]['parts'][a]['mus'] if a != 'all' else self.pieces[x].mu()

    def sub(self, x, a, y, b):
        key = (x, a, y, b)
        if key not in self._sub:
            A = hom_system(self.ph(x, a), self.ph(y, b), self.pmus(x, a), self.pmus(y, b), self.p)
            self._sub[key] = A.shape[1] - rank_mod(A, self.p)
        return self._sub[key]


def run_family(name, d, pidx=None):
    t0 = time.time()
    n, w = parse_family(name)
    S = family_S(d, n, w)
    L = lcm(d, n)
    U = units(L)
    primes = pick_primes(d, n)
    if pidx is not None:
        primes = [primes[pidx]]
    log = []
    PD = [PrimeData(name, d, p) for p in primes]
    P0 = PD[0].pieces
    keys = sorted(P0)
    h = {x: P0[x].h for x in keys}
    for pd in PD[1:]:
        assert sorted(pd.pieces) == keys and all(pd.pieces[x].h == h[x] for x in keys)
    for x in keys:
        assert h[x] == tuple(PD[0].R.koszul(e, P0[x].cc) if e >= 0 else 0 for e in P0[x].deg)
    dim = {x: sum(h[x]) for x in keys}
    pure = {x: h[x][0] == 0 and h[x][2] == 0 for x in keys}
    orb = {x: frozenset(gal(t, x, d, n) for t in U) for x in keys}
    for x in keys:
        assert all(y in P0 and dim[y] == dim[x] for y in orb[x])
    pure_orbit = {x: all(pure[y] for y in orb[x]) for x in keys}
    conj = {x: ((-x[0]) % d, (-x[1]) % n) for x in keys}
    planes = curve_planes(d, n, w, S)
    delta = {x: curve_delta(planes, x[0], x[1], d, n) for x in keys}
    curve_np = {x: delta[x] > 0 and not pure_orbit[x] for x in keys}
    vps = valid_perms(d, n, w, S)
    valid = [v for v in vps if v['ii'] and v['iii']]
    invol = [v for v in valid if v['order'] == 2]
    gm = gmax_check(d, n, w, S)
    usw = u_swap_flags(d, S)
    # automorphism orbits (Gamma acts trivially on labels; valid perms permute labels)
    Aorb = {}
    for x in keys:
        seen = {x}
        stack = [x]
        while stack:
            z = stack.pop()
            for v in valid:
                y = (z[0], v['cmap'][z])
                if y not in seen:
                    seen.add(y); stack.append(y)
        Aorb[x] = seen
    # ---------------- certificates
    comm = {x: [pd.comm(x) for pd in PD] for x in keys}
    certE, how = {}, {}
    for x in keys:
        certE[x] = min(comm[x]) == 1
        how[x] = 'comm=1' if certE[x] else ''
    hull = {}
    def hull_of(x):
        if x not in hull:
            hull[x] = max(pd.hull(x) for pd in PD) if h[x][0] else 0
        return hull[x]
    def l4p_any(x, h11eff):
        for k, pd in enumerate(PD):
            if l4p_test_eff(pd.pieces[x], random.Random(f'Z3l4p|{name}|{d}|{x}|{k}'), h11eff):
                return True
        return False
    for x in keys:
        if certE[x] or pure[x] or delta[x] > 0:
            continue
        c = conj[x]
        if h[x][0] == 1 and hull_of(x) == dim[x]:
            certE[x], how[x] = True, 'hull,h20=1 (L2)'
        elif h[x][2] == 1 and hull_of(c) == dim[x]:
            certE[x], how[x] = True, 'conj hull,h02=1 (L2+0.5)'
        elif h[x][0] == 2 and hull_of(x) == dim[x] and l4p_any(x, h[x][1]):
            certE[x], how[x] = True, "hull+L4'"
        elif h[x][2] == 2 and hull_of(c) == dim[x] and l4p_any(c, h[c][1]):
            certE[x], how[x] = True, "conj hull+L4'"
    certU, hU, errors = {}, {}, []
    for x in keys:
        certU[x] = False
        if not curve_np[x]:
            continue
        dl, c = delta[x], conj[x]
        hU[x] = (h[x][0], h[x][1] - dl, h[x][2])
        if h[x][0] >= 1:
            hl = hull_of(x)
            if hl > dim[x] - dl:
                errors.append(f'hull {hl} > dim-delta at {x}')
        if h[x][0] == 1 and hull_of(x) == dim[x] - dl:
            certU[x], how[x] = True, 'U: hull=dim-delta,h20=1 (P1.7a)'
        elif h[x][2] == 1 and hull_of(c) == dim[x] - dl:
            certU[x], how[x] = True, 'U: conj hull=dim-delta,h02=1 (P1.7b)'
        elif h[x][0] == 2 and hull_of(x) == dim[x] - dl and l4p_any(x, h[x][1] - dl):
            certU[x], how[x] = True, "U: hull+L4' (P1.7c)"
        elif h[x][2] == 2 and hull_of(c) == dim[x] - dl and l4p_any(c, h[c][1] - dl):
            certU[x], how[x] = True, "U: conj hull+L4' (P1.7c)"
    # rigidity of orbits by Hodge numbers (S2 L3(iii)) with an irreducible member
    def rigid(x):
        O = orb[x]
        if not any(certE[y] and not pure[y] for y in O):
            return False
        for t in U:
            if gal(t, x, d, n) == x:
                continue
            if not any(h[gal(t, s, d, n)] != h[s] for s in O):
                return False
        return True
    rig = {x: rigid(x) for x in keys}
    # ---------------- split data
    SD = {}
    for v in invol:
        SD[v['pi']] = [SplitData(name, d, p, v, keys) for p in primes]
    def stable_under(x, pi):
        return x in SD[pi][0].stable
    part_info = {}
    for pi, sds in SD.items():
        sd0 = sds[0]
        for x in sd0.stable:
            for a in (1, -1):
                hp = sd0.ph(x, a)
                part_info[(pi, x, a)] = dict(h=hp, pure=(hp[0] == 0 and hp[2] == 0), dim=sum(hp))
        for x in sd0.stable:
            for a in (1, -1):
                pinf = part_info[(pi, x, a)]
                pinf['pure_orbit'] = all(part_info[(pi, y, a)]['pure'] for y in orb[x])
                sc = [sd.sub(x, a, x, a) for sd in sds] if pinf['dim'] else [0, 0]
                pinf['comm'] = sc
                pinf['cert'] = pinf['dim'] > 0 and min(sc) == 1
    def pinfo(pi, x, a):
        if a == 'all':
            return dict(h=h[x], pure=pure[x], dim=dim[x], pure_orbit=pure_orbit[x], cert=certE[x])
        return part_info[(pi, x, a)]
    def split_lb(x, y, pi):
        sd = SD[pi][0]
        tot = 0
        for a in sd.parts(x):
            for b in sd.parts(y):
                A, B = pinfo(pi, x, a), pinfo(pi, y, b)
                if A['dim'] == 0 or B['dim'] == 0:
                    continue
                if A['pure_orbit'] and B['pure_orbit']:
                    tot += A['dim'] * B['dim']
                elif x == y and a == b:
                    tot += 1
        return tot
    def split_ub(x, y, pi, k):
        sd = SD[pi][k]
        tot = 0
        for a in sd.parts(x):
            for b in sd.parts(y):
                A, B = pinfo(pi, x, a), pinfo(pi, y, b)
                if A['dim'] == 0 or B['dim'] == 0:
                    continue
                if A['pure_orbit'] and B['pure_orbit']:
                    u = A['dim'] * B['dim']
                elif A['pure'] and B['pure']:
                    u = A['dim'] * B['dim']
                elif (A['cert'] and not A['pure'] and B['pure']) or (B['cert'] and not B['pure'] and A['pure']):
                    u = 0
                elif A['cert'] and B['cert'] and A['h'] != B['h']:
                    u = 0
                elif x == y and a == b and A['cert']:
                    u = 1
                else:
                    u = sd.sub(x, a, y, b)
                tot += u
        return tot
    # ---------------- lower bounds
    def ns_dim(x):
        if pure_orbit[x]:
            return dim[x]
        if curve_np[x]:
            return delta[x]
        return 0
    def lb_member(x, y):
        nx, ny = ns_dim(x), ns_dim(y)
        g = 1 if (y in Aorb[x] and dim[x] > nx) else 0
        best, src = nx * ny + g, 'Gamma/A graphs + NSxNS'
        for pi in SD:
            if stable_under(x, pi) or stable_under(y, pi):
                v = split_lb(x, y, pi)
                if v > best:
                    best, src = v, 'split graphs + NSxNS'
        return best, src
    # ---------------- structural upper bounds at a member
    def member_ub(x, y):
        c = []
        if pure_orbit[x] and pure_orbit[y]:
            return dim[x] * dim[y], 'NSxNS'
        if pure_orbit[y] and certU[x]: c.append((delta[x] * dim[y], 'curve-b(P1.8b)'))
        if pure_orbit[x] and certU[y]: c.append((delta[y] * dim[x], 'curve-b(P1.8b)'))
        if certE[x] and not pure[x] and pure[y]: c.append((0, 'C5.2v'))
        if certE[y] and not pure[y] and pure[x]: c.append((0, 'C5.2v'))
        if certU[x] and pure[y]: c.append((delta[x] * dim[y], 'curve-v'))
        if certU[y] and pure[x]: c.append((delta[y] * dim[x], 'curve-v'))
        if x == y:
            if certE[x]: c.append((1, 'Schur'))
            if certU[x]: c.append((1 + delta[x] ** 2, 'curve-a(P1.8a)'))
        else:
            if certE[x] and certE[y]:
                c.append((0, 'C5.2iv') if h[x] != h[y] else (1, 'Schur-iso'))
            if certE[x] and any(h[x][i] > h[y][i] for i in range(3)): c.append((0, 'C5.2i(inj)'))
            if certE[y] and any(h[y][i] > h[x][i] for i in range(3)): c.append((0, 'C5.2ii(surj)'))
            # P-component killed by C5.2(v): needs the certified E-side to be NON-pure (fix after self-review)
            if certU[x] and certE[y] and not pure[y]: c.append((0 if hU[x] != h[y] else 1, 'curve-c(U irr)'))
            if certE[x] and not pure[x] and certU[y]: c.append((0 if h[x] != hU[y] else 1, 'curve-d(U irr)'))
            if certU[x] and not certE[y] and not certU[y] and not pure[y]:
                if any(hU[x][i] > h[y][i] for i in range(3)): c.append((delta[x] * dim[y], 'curve-inj'))
            if certU[y] and not certE[x] and not certU[x] and not pure[x]:
                if any(hU[y][i] > h[x][i] for i in range(3)): c.append((delta[y] * dim[x], 'curve-surj'))
        return min(c) if c else (INF, None)
    # ---------------- block orbits
    pairs = [(x, y) for x in keys for y in keys if x[0] == y[0]]
    seen = set()
    borbs = []
    for (x, y) in pairs:
        if (x, y) in seen:
            continue
        mem = sorted(set((gal(t, x, d, n), gal(t, y, d, n)) for t in U))
        seen.update(mem)
        borbs.append(mem)
    blocks = []
    for mem in borbs:
        x0, y0 = mem[0]
        cands = list(mem) + [(y, x) for (x, y) in mem if (y, x) not in mem]
        rec = dict(rep=[list(x0), list(y0)], size=len(mem), h=[h[x0], h[y0]], details=[])
        lb, lbsrc = max(lb_member(x, y) for (x, y) in cands)
        ub, tool = INF, None
        for (x, y) in cands:
            u, t = member_ub(x, y)
            if u < ub:
                ub, tool = u, f'{t}@{x}->{y}'
        if ub > lb:
            if x0 != y0 and orb[x0] == orb[y0] and rig[x0]:
                ub, tool = 0, 'R3b(L3 rigid orbit)'
        if ub > lb and ub > 1:
            if any(certE[z] for z in orb[x0]) and any(certE[z] for z in orb[y0]):
                ub, tool = 1, 'cap1(P4.2e)'
        if ub > lb:
            sel = sorted(cands, key=lambda q: (pure[q[0]] + pure[q[1]], dim[q[0]] * dim[q[1]]))[:MAX_HOM]
            for (x, y) in sel:
                for k, pd in enumerate(PD):
                    v = pd.hom(x, y)
                    rec['details'].append(f'homcomm{x}->{y}@p{k}={v}')
                    if v < ub:
                        ub, tool = v, f'homcomm@{x}->{y}'
                if ub <= lb:
                    break
        if ub > lb and SD:
            sel = [q for q in cands if any(stable_under(q[0], pi) or stable_under(q[1], pi) for pi in SD)][:MAX_SPLIT]
            for (x, y) in sel:
                for pi in SD:
                    if not (stable_under(x, pi) or stable_under(y, pi)):
                        continue
                    for k in range(len(primes)):
                        v = split_ub(x, y, pi, k)
                        rec['details'].append(f'split{pi}{x}->{y}@p{k}={v}')
                        if v < ub:
                            ub, tool = v, f'split(T3)@{x}->{y}'
                if ub <= lb:
                    break
        if ub > lb:
            sel = [q for q in cands if h[q[0]][0] >= 1]
            sel.sort(key=lambda q: (not (delta[q[0]] > 0 and certE[q[1]] and not pure[q[1]]), pure[q[1]], dim[q[0]] * dim[q[1]]))
            for (x, y) in sel[:MAX_J]:
                for k, pd in enumerate(PD):
                    J = pd.J(x, y)
                    if delta[x] > 0 and certE[y] and not pure[y]:
                        v, t = J - delta[x] * dim[y], 'J_V(L1.9)'
                    else:
                        v, t = J, 'J(T2)'
                    rec['details'].append(f'J{x}->{y}@p{k}={J}' + (f'(J_V={v})' if t.startswith('J_V') else ''))
                    if v < ub:
                        ub, tool = v, f'{t}@{x}->{y}'
                if ub <= lb:
                    break
        if lb > ub:
            errors.append(f'lb {lb} > ub {ub} at block {mem[0]}')
        rec.update(lb=lb, lbsrc=lbsrc, ub=ub, tool=tool, closed=(lb == ub))
        blocks.append(rec)
    Ylb = sum(b['size'] * b['lb'] for b in blocks)
    Yub = sum(b['size'] * min(b['ub'], 10 ** 6) for b in blocks)
    res = dict(name=name, d=d, n=n, w=list(w), S=len(S), primes=primes, redraws=[pd.tries for pd in PD],
               gmax=list(gm), uswap=usw, planes=[dict(i=pl['i'], M=list(pl['M']), g=pl['g'], e=pl['e']) for pl in planes],
               perms=[dict(pi=list(v['pi']), ii=v['ii'], iii=v['iii'], order=v['order']) for v in vps],
               split_checks={str(pi): [[list(c[0]), c[1], c[2]] for sd in sds for c in sd.checks] for pi, sds in SD.items()},
               split_redraws={str(pi): [sd.tries for sd in sds] for pi, sds in SD.items()},
               pieces={str(x): dict(h=h[x], dim=dim[x], pure=pure[x], pure_orbit=pure_orbit[x], orbit=len(orb[x]),
                                    delta=delta[x], comm=comm[x], certE=certE[x], certU=certU[x], how=how[x],
                                    hull=hull.get(x), rigid=rig[x],
                                    split={str(pi): {str(a): part_info[(pi, x, a)] for a in (1, -1)}
                                           for pi in SD if (pi, x, 1) in part_info}) for x in keys},
               blocks=blocks, Ylb=Ylb, Yub=Yub, closed=all(b['closed'] for b in blocks), errors=errors,
               seconds=round(time.time() - t0, 1))
    return res


def l4p_test_eff(P, rng, h11eff, trials=4):
    """L4' test with the threshold computed for h11eff (= h11(U) for a curve-split piece)."""
    R = P.R
    p = R.p
    if P.h[0] != 2:
        return False
    k1 = h11eff // 2 + 1
    S = R.S
    if k1 > min(P.h[1], len(S)):
        return False
    mus = P.mu()
    M1 = np.column_stack([mus[i][0][:, 0] for i in range(len(S))]) % p
    M2 = np.column_stack([mus[i][0][:, 1] for i in range(len(S))]) % p
    g, lead = None, False
    for _ in range(trials):
        A = np.array([[rng.randrange(p) for _ in range(P.h[1])] for _ in range(k1)], dtype=np.int64)
        B = np.array([[rng.randrange(p) for _ in range(k1)] for _ in range(len(S))], dtype=np.int64)
        X1 = ((A @ M1) % p @ B) % p
        X2 = ((A @ M2) % p @ B) % p
        xs = list(range(k1 + 1))
        ys = [det_mod((s * X1 + X2) % p, p) for s in xs]
        f = nmod_poly(interp_mod(xs, ys, p), p)
        if det_mod(X1, p) % p:
            lead = True
        g = f if g is None else g.gcd(f)
    return g is not None and g.degree() == 0 and lead


def fmt(res):
    lines = []
    lines.append(f"== {res['name']} d={res['d']} n={res['n']} w={res['w']} |S|={res['S']} primes={res['primes']} "
                 f"redraws={res['redraws']} |H|,|Gamma img|={res['gmax']} u-swap idx={res['uswap']} ({res['seconds']}s)")
    lines.append(f"   planes={res['planes']}")
    lines.append(f"   perms={res['perms']} split_redraws={res['split_redraws']}")
    lines.append(f"   split checks (Pi^2=1, Pi mu=mu Pi) all true: "
                 f"{all(c[1] and c[2] for v in res['split_checks'].values() for c in v)}")
    lines.append(f"   Y lb={res['Ylb']} ub={res['Yub']} closed={res['closed']} errors={res['errors']}")
    lines.append('   pieces (a,chi): h dim pure/pureorbit orbit delta comm cert how hull rigid split')
    for k, v in res['pieces'].items():
        sp = {pi: {a: (q['h'], q['pure_orbit'], q['comm'], q['cert']) for a, q in d_.items()} for pi, d_ in v['split'].items()}
        lines.append(f"     {k}: h={tuple(v['h'])} dim={v['dim']} pure={int(v['pure'])}/{int(v['pure_orbit'])} orb={v['orbit']} "
                     f"delta={v['delta']} comm={v['comm']} E={int(v['certE'])} U={int(v['certU'])} [{v['how']}] hull={v['hull']} "
                     f"rigid={int(v['rigid'])} {sp if sp else ''}")
    lines.append('   blocks: rep size h lb ub closed tool | lbsrc | details')
    for b in res['blocks']:
        lines.append(f"     {tuple(b['rep'][0])}->{tuple(b['rep'][1])} x{b['size']} h={tuple(b['h'][0])}->{tuple(b['h'][1])} "
                     f"lb={b['lb']} ub={b['ub']} {'CLOSED' if b['closed'] else 'OPEN'} {b['tool']} | {b['lbsrc']} | "
                     + ' '.join(b['details']))
    return '\n'.join(lines)


def work(fd):
    name, d, pidx = fd
    try:
        res = run_family(name, d, pidx)
        tag = f"{name}_d{d}" + ('' if pidx is None else f'_only_p{pidx}')
        with open(os.path.join(OUTDIR, f'z3_{tag}.json'), 'w', encoding='utf-8') as f:
            json.dump(res, f, default=str)
        with open(os.path.join(OUTDIR, f'z3_{tag}.txt'), 'w', encoding='utf-8') as f:
            f.write(fmt(res) + '\n')
        return (name, d, res['Ylb'], res['Yub'], res['closed'], res['errors'], res['seconds'])
    except Exception:
        return (name, d, None, None, False, [traceback.format_exc()], 0)


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    os.makedirs(OUTDIR, exist_ok=True)
    args = sys.argv[1:]
    pidx = None
    if args and args[0].startswith('--p'):
        pidx = int(args[0][3:])
        args = args[1:]
    nproc = int(args[0]) if args else 1
    if len(args) > 1:
        todo = []
        for s in args[1:]:
            a, b = s.split(':')
            todo.append((a, int(b), pidx))
        tag = 'partial'
    else:
        todo = [(a, b, pidx) for (a, b) in FAMILIES + CONTROLS]
        tag = 'all' if pidx is None else f'only_p{pidx}'
    todo.sort(key=lambda fd: -(fd[1] ** 3 * parse_family(fd[0])[0]))
    if nproc > 1:
        import multiprocessing as mp
        with mp.Pool(nproc) as pool:
            out = pool.map(work, todo, chunksize=1)
    else:
        out = [work(fd) for fd in todo]
    lines = []
    for r in out:
        lines.append(f"{r[0]:14s} d={r[1]:2d}  Y lb={r[2]} ub={r[3]} closed={r[4]} errors={r[5]} {r[6]}s")
    txt = '\n'.join(lines)
    print(txt)
    with open(os.path.join(OUTDIR, f'z3_summary_{tag}.txt'), 'w', encoding='utf-8') as f:
        f.write(txt + '\n')
