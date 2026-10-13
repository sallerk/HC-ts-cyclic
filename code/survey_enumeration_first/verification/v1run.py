# V1 family runner: pieces, Galois orbits, automorphisms, rigidity certificates, block classification,
# Y (Hodge classes in sum_j End(H_j)) and C (rank End J(C_F)) for the very general member.
# Independent implementation from S2_lemmas.md (Theorem T, L2, L3, L4', L5, L6) and the rule list S1_tests.md.
import os, sys, json, time, random, itertools, traceback
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v1lib import *
from math import gcd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTD = os.path.join(BASE, 'outputs')


def family_monomials(d, w, n):
    return [m for m in monos(d) if (w[0] * m[0] + w[1] * m[1] + w[2] * m[2]) % n == 0]


def kstr(key):
    return ':'.join(str(x) for x in key)


class Fam:
    def __init__(self, name, d, w, n, S=None, nprimes=2, start=50000):
        self.name, self.d, self.w, self.n = name, d, tuple(w), n
        self.S = list(S) if S else family_monomials(d, w, n)
        self.sw = sum(w) % n
        self.perms = [pi for pi in allowed_perms(self.S) if pi != (0, 1, 2)]
        self.gmax = gmax_order(self.S)
        self.primes = pick_primes(lcm(8, d, n), nprimes, start=start)
        self.rings = {}
        for idx, p in enumerate(self.primes):
            rng = random.Random(f'V1|{name}|d{d}|{p}|{idx}')
            F = {m: rng.randrange(1, p) for m in self.S}
            R = Ring(d, w, n, p, F)
            bad = R.check_hilbert()
            if bad:
                raise RuntimeError(f'Hilbert/smoothness check failed at p={p}: {bad[:5]}')
            self.rings[p] = R
        self.Tk = [{(0, 0, 0)}]
        for k in range(3):
            self.Tk.append({madd(a, b) for a in self.Tk[-1] for b in self.S})
        self.pieces = {}
        R0 = self.rings[self.primes[0]]
        for a in range(1, d):
            shift = 3 + a
            for c in range(n):
                h = tuple(R0.dimR(q * d - shift, c) for q in (1, 2, 3))
                for p in self.primes[1:]:
                    assert h == tuple(self.rings[p].dimR(q * d - shift, c) for q in (1, 2, 3))
                if sum(h):
                    self.pieces[('S', a, (c + self.sw) % n)] = dict(shift=shift, c=c, h=h, qmax=3)
        for c in range(n):
            h = (R0.dimR(d - 3, c), R0.dimR(2 * d - 3, c))
            if sum(h):
                self.pieces[('C', (c + self.sw) % n)] = dict(shift=3, c=c, h=h, qmax=2)
        self.NS = lcm(d, n)
        self.unitsS = [k for k in range(1, self.NS) if gcd(k, self.NS) == 1]
        self.unitsC = [k for k in range(1, n) if gcd(k, n) == 1] or [1]
        # Hodge symmetry under complex conjugation (tests the character convention incl. Omega)
        self.sym_ok = True
        for key, P in self.pieces.items():
            ck = self.act(-1, key)
            if ck not in self.pieces or self.pieces[ck]['h'] != tuple(reversed(P['h'])):
                self.sym_ok = False
        # automorphism action on characters
        self.permmaps = []
        self.perm_ok = True
        for pi in self.perms:
            mp = {}
            for key, P in self.pieces.items():
                ws = set()
                for q in range(1, P['qmax'] + 1):
                    for m in R0.cols(q * d - P['shift'], P['c'])[0]:
                        ws.add(R0.wt(perm_apply(pi, m)))
                if len(ws) != 1:
                    self.perm_ok = False
                    continue
                c2 = ws.pop()
                k2 = (key[0], key[1], (c2 + self.sw) % n) if key[0] == 'S' else ('C', (c2 + self.sw) % n)
                if k2 not in self.pieces or self.pieces[k2]['h'] != P['h']:
                    self.perm_ok = False
                mp[key] = k2
            self.permmaps.append(mp)
        # K0 = span of x_i d_j F with w_i = w_j (exact kernel of every theta)
        self.K0 = {}
        idx = {m: k for k, m in enumerate(self.S)}
        for p, R in self.rings.items():
            vecs = []
            for i in range(3):
                for j in range(3):
                    if (self.w[i] - self.w[j]) % n == 0:
                        v = [0] * len(self.S)
                        for m, cf in R.dF[j].items():
                            mm = list(m)
                            mm[i] += 1
                            v[idx[tuple(mm)]] = (v[idx[tuple(mm)]] + cf) % p
                        vecs.append(v)
            self.K0[p] = rank_rows(vecs, len(self.S), p)
        self._theta, self._hull, self._kern, self._ub1, self._l4 = {}, {}, {}, {}, {}
        self.rng = random.Random(f'V1-cert|{name}|{d}')

    # ---------------- group actions
    def act(self, k, key):
        if key[0] == 'S':
            return ('S', (k * key[1]) % self.d, (k * key[2]) % self.n)
        return ('C', (k * key[1]) % self.n)

    def units(self, T):
        return self.unitsS if T == 'S' else self.unitsC

    def gorbit(self, key):
        return tuple(sorted(set(self.act(k, key) for k in self.units(key[0]))))

    def aorbit(self, key):
        seen, todo = {key}, [key]
        while todo:
            x = todo.pop()
            for mp in self.permmaps:
                y = mp.get(x)
                if y is not None and y not in seen:
                    seen.add(y)
                    todo.append(y)
        return seen

    def lb_aut(self, B):
        """lower bound for a block containing graph classes: max over (x,y) in B, y in A.x, of the largest orbit of
        Stab_P(x) on the monomial basis of the top Hodge piece of E_x (generalised permutation matrices with pairwise
        different images of one basis vector are linearly independent)."""
        best = 0
        for x, y in B:
            if y not in self.aorbit(x):
                continue
            best = max(best, 1)
            P = self.pieces[x]
            R0 = self.rings[self.primes[0]]
            top = R0.cols(self.d - P['shift'], P['c'])[0]
            stab = [None] + [pi for pi, mp in zip(self.perms, self.permmaps) if mp.get(x) == x]
            for m in top:
                orbit = {m} | {perm_apply(pi, m) for pi in stab if pi is not None}
                best = max(best, len(orbit))
        return best

    def h(self, key):
        return self.pieces[key]['h']

    # ---------------- per-prime data
    def theta(self, p, key):
        if (p, key) not in self._theta:
            P = self.pieces[key]
            self._theta[(p, key)] = self.rings[p].theta(P['shift'], P['c'], self.S, P['qmax'])
        return self._theta[(p, key)]

    def hull(self, p, key):
        if (p, key) not in self._hull:
            P = self.pieces[key]
            self._hull[(p, key)] = self.rings[p].hull(P['shift'], P['c'], self.S, self.Tk, P['qmax'])[0]
        return self._hull[(p, key)]

    def l4(self, p, key):
        if (p, key) not in self._l4:
            h = self.h(key)
            th = self.theta(p, key)[0]
            M1 = [[th[vi][r][0] for vi in range(len(self.S))] for r in range(h[1])]
            M2 = [[th[vi][r][1] for vi in range(len(self.S))] for r in range(h[1])]
            k = h[1] // 2 + 1
            cert = minrank_cert(M1, M2, k, p, self.rng)
            lowb = minrank_lower(M1, M2, p, self.rng)
            self._l4[(p, key)] = dict(k=k, cert=cert, minrank_ge=lowb)
        return self._l4[(p, key)]

    def kern(self, p, key, l):
        if (p, key, l) not in self._kern:
            th = self.theta(p, key)[l]
            rows = [[x for r in M for x in r] for M in th]
            ncol = len(rows[0]) if rows else 0
            K = left_kernel(rows, ncol, p)
            rk = len(self.S) - len(K)
            valid = rk == min(len(self.S) - self.K0[p], ncol)
            self._kern[(p, key, l)] = (rk, tuple(K), valid)
        return self._kern[(p, key, l)]

    def ub1(self, p, s, t):
        if (p, s, t) not in self._ub1:
            self._ub1[(p, s, t)] = ub1_block(self.theta(p, s), self.theta(p, t), self.h(s), self.h(t), p)
        return self._ub1[(p, s, t)]

    def ub1min(self, B, lb):
        cost = lambda st: sum(a * b for a, b in zip(self.h(st[0]), self.h(st[1])))
        best = None
        vals = []
        for st in sorted(B, key=cost):
            for p in self.primes:
                v = self.ub1(p, st[0], st[1])
                vals.append((kstr(st[0]), kstr(st[1]), p, v))
                best = v if best is None else min(best, v)
                if best <= lb:
                    return best, vals
        return best, vals

    # ---------------- orbit kinds
    def classify_orbits(self, T):
        pure = 'pure' if T == 'S' else 'cm'
        keys = sorted(k for k in self.pieces if k[0] == T)
        orbits, self_orb = [], {}
        for key in keys:
            if key in self_orb:
                continue
            O = self.gorbit(key)
            for x in O:
                assert x in self.pieces and sum(self.h(x)) == sum(self.h(key))
                self_orb[x] = O
            orbits.append(O)
        info = {}
        for O in orbits:
            if T == 'S':
                ispure = all(self.h(x)[0] == 0 and self.h(x)[2] == 0 for x in O)
            else:
                ispure = all(self.h(x)[0] * self.h(x)[1] == 0 for x in O)
            rec = dict(members=[kstr(x) for x in O], h={kstr(x): self.h(x) for x in O})
            if ispure:
                rec['kind'] = pure
            else:
                cands = sorted([x for x in O if self.h(x)[0] in (1, 2)], key=lambda x: (self.h(x)[0], x))
                cert, attempts = None, []
                for x in cands:
                    for p in self.primes:
                        hr = self.hull(p, x)
                        full = hr[-1] == sum(self.h(x))
                        a = dict(member=kstr(x), p=p, hull=hr, full=full)
                        ok = full
                        if full and self.h(x)[0] == 2:
                            a['l4'] = self.l4(p, x)
                            ok = a['l4']['cert']
                        attempts.append(a)
                        if ok:
                            cert = x
                            break
                    if cert:
                        for p in self.primes:   # report the other prime as well
                            if not any(a['p'] == p and a['member'] == kstr(x) for a in attempts):
                                a = dict(member=kstr(x), p=p, hull=self.hull(p, x))
                                if self.h(x)[0] == 2:
                                    a['l4'] = self.l4(p, x)
                                attempts.append(a)
                        break
                rec['attempts'] = attempts
                rec['kind'] = 'rigid' if cert else 'res'
                rec['cert'] = cert
            info[O] = rec
        # Galois rigidity (End_Hdg(H_O) = K, L3(iii)) for rigid orbits
        for O, rec in info.items():
            if rec['kind'] != 'rigid':
                continue
            s0 = rec['cert']
            gr, why = True, []
            for t in sorted(set(self.act(k, s0) for k in self.units(T))):
                if t == s0:
                    continue
                if self.h(t) != self.h(s0):
                    why.append((kstr(t), 'hodge'))
                    continue
                if t in self.aorbit(s0):
                    gr = False
                    why.append((kstr(t), 'automorphism'))
                    break
                if self.r5valid(s0, t):
                    why.append((kstr(t), 'higgs-kernel'))
                    continue
                Bp = sorted(set((self.act(k, s0), self.act(k, t)) for k in self.units(T)))
                if self.rule2(Bp, 'rigid', 'rigid'):
                    why.append((kstr(t), 'hodge-orbit-sum'))
                    continue
                v, _ = self.ub1min(Bp, 0)
                if v == 0:
                    why.append((kstr(t), 'ub1=0'))
                    continue
                gr = False
                why.append((kstr(t), 'undecided'))
                break
            rec['galois_rigid'] = gr
            rec['grigid_why'] = why
        return orbits, self_orb, info

    def r5valid(self, s, t):
        if self.h(s) != self.h(t):
            return False
        for p in self.primes:
            for l in range(self.pieces[s]['qmax'] - 1):
                rs, Ks, vs = self.kern(p, s, l)
                rt, Kt, vt = self.kern(p, t, l)
                if vs and vt and Ks != Kt:
                    return True
        return False

    def r5naive(self, s, t):
        if self.h(s) != self.h(t):
            return False
        for l in range(self.pieces[s]['qmax'] - 1):
            if all(self.kern(p, s, l)[1] != self.kern(p, t, l)[1] for p in self.primes):
                return True
        return False

    # ---------------- blocks
    def blocks(self, T):
        orbits, orb, info = self.classify_orbits(T)
        pure = 'pure' if T == 'S' else 'cm'
        keys = sorted(k for k in self.pieces if k[0] == T)
        kind = lambda x: info[orb[x]]['kind']
        done, recs = set(), []
        conflicts = []
        for s in keys:
            for t in keys:
                if T == 'S' and s[1] != t[1]:
                    continue
                if (s, t) in done:
                    continue
                B = sorted(set((self.act(k, s), self.act(k, t)) for k in self.units(T)))
                done.update(B)
                hs, ht = self.h(s), self.h(t)
                ds, dt = sum(hs), sum(ht)
                Ks, Kt = kind(s), kind(t)
                aut = t in self.aorbit(s)
                rec = dict(s=kstr(s), t=kstr(t), size=len(B), hs=hs, ht=ht, ks=Ks, kt=Kt, aut=aut)
                if Ks == pure and Kt == pure:
                    if T == 'S':
                        lb = ub = ds * dt
                        rule = 'R1:NSxNS'
                    else:
                        typ = lambda x: 0 if self.h(x)[1] == 0 else 1
                        matched = all(typ(x) == typ(y) for x, y in B)
                        lb = ub = ds * dt if matched else 0
                        rule = 'R1:CMxCM'
                elif (Ks == pure and Kt == 'rigid') or (Ks == 'rigid' and Kt == pure):
                    lb = ub = 0
                    rule = 'R1:pure-rigid'
                    if aut:
                        conflicts.append(rec)
                elif Ks == pure or Kt == pure:
                    rule = 'R1:pure-res'
                    lb = self.lb_aut(B) if aut else 0
                    ub, vals = self.ub1min(B, lb)
                    rec['ub1'] = vals
                else:
                    zr = []
                    if self.rule2(B, Ks, Kt):
                        zr.append('R2')
                    if self.rule2c(s, t, Ks, Kt, info, orb):
                        zr.append('R2c')
                    both = Ks == 'rigid' and Kt == 'rigid'
                    if both and info[orb[s]]['galois_rigid'] and info[orb[t]]['galois_rigid']:
                        if orb[s] != orb[t] and not (len(B) == len(orb[s]) == len(orb[t])):
                            zr.append('R3')
                        if orb[s] == orb[t] and s != t:
                            zr.append('R3b')
                    # OLD rule 5: exact (rank-certified) kernels compared at the pair (s,t) itself, both orbits rigid.
                    # FIXED rule 5' (referee P1/P2): only at transported pairs with a CERTIFIED embedding:
                    # certified source x: K(E_y) in K(E_x); certified target y: K(E_x) in K(E_y); violation by an
                    # exact (rank-certified) kernel => block 0.
                    r5old = both and not aut and self.r5valid(s, t)
                    r5new = (not aut) and self.r5onesided(B, Ks, Kt, info, orb)
                    if both and not aut:
                        rec['r5naive'] = self.r5naive(s, t)
                    rec['r5old'], rec['r5new'] = r5old, r5new
                    zo = zr + (['R5'] if r5old else [])
                    zf = zr + (["R5'"] if r5new else [])
                    rec['zero_rules'], rec['zero_rules_fixed'] = zo, zf

                    def resid(store=True):
                        u, vals = self.ub1min(B, 0)
                        rec['ub1' if store else 'ub1_check'] = vals
                        return min(u, 1) if both else u
                    if aut:
                        lb = self.lb_aut(B)
                        rec['lb_aut'] = lb
                        if zo or zf or (both and lb > 1):
                            conflicts.append(rec)
                        if both:
                            ub, rule = 1, 'R4:aut'
                        else:
                            ub, vals = self.ub1min(B, lb)
                            rec['ub1'] = vals
                            rule = 'R6:aut-res'
                        fixed = (rule, lb, ub)
                    else:
                        lb = 0
                        if zo:
                            ub, rule = 0, zo[0]
                        else:
                            ub, rule = resid(), 'R6:res'
                        if zf:
                            fixed = (zf[0], 0, 0)
                        else:
                            fixed = ('R6:res', 0, ub if not zo else resid(False))
                        if zo and set(zo) <= {'R5', 'R3', 'R3b'}:
                            rec['only_r5_or_grigid'] = True
                            rec['ub1_if_rule_dropped'] = fixed[2] if not zf else resid(False)
                    rec['fixed'] = fixed
                rec.update(rule=rule, lb=lb, ub=ub, closed=(lb == ub))
                if 'fixed' not in rec:
                    rec['fixed'] = (rule, lb, ub)
                rec['closed_fixed'] = rec['fixed'][1] == rec['fixed'][2]
                recs.append(rec)
        return info, recs, conflicts

    def kincl_viol(self, a, b):
        """True if some exact (rank-certified) kernel K^l(E_a) is not contained in the mod-p kernel K^l(E_b)."""
        for p in self.primes:
            for l in range(self.pieces[a]['qmax'] - 1):
                ra, Ka, va = self.kern(p, a, l)
                rb, Kb, vb = self.kern(p, b, l)
                if not va or not Ka:
                    continue
                if rank_rows(list(Kb) + list(Ka), len(self.S), p) > len(Kb):
                    return True
        return False

    def r5onesided(self, B, Ks, Kt, info, orb):
        s0 = info[orb[B[0][0]]]['cert'] if Ks == 'rigid' else None
        t0 = info[orb[B[0][1]]]['cert'] if Kt == 'rigid' else None
        for x, y in B:
            if s0 is not None and x == s0 and self.kincl_viol(y, x):
                return True
            if t0 is not None and y == t0 and self.kincl_viol(x, y):
                return True
        return False

    def rule2(self, B, Ks, Kt):
        nl = len(self.h(B[0][0]))
        if Ks == 'rigid':
            for x in set(a for a, _ in B):
                ys = [y for a, y in B if a == x]
                if any(self.h(x)[l] > sum(self.h(y)[l] for y in ys) for l in range(nl)):
                    return True
        if Kt == 'rigid':
            for y in set(b for _, b in B):
                xs = [x for x, b in B if b == y]
                if any(self.h(y)[l] > sum(self.h(x)[l] for x in xs) for l in range(nl)):
                    return True
        return False

    def rule2c(self, s, t, Ks, Kt, info, orb):
        T = s[0]
        nl = len(self.h(s))
        if Ks == 'rigid':
            s0 = info[orb[s]]['cert']
            for k in self.units(T):
                if self.act(k, s) == s0:
                    t2 = self.act(k, t)
                    if any(self.h(s0)[l] > self.h(t2)[l] for l in range(nl)):
                        return True
        if Kt == 'rigid':
            t0 = info[orb[t]]['cert']
            for k in self.units(T):
                if self.act(k, t) == t0:
                    s2 = self.act(k, s)
                    if any(self.h(s2)[l] < self.h(t0)[l] for l in range(nl)):
                        return True
        if Ks == 'rigid' and Kt == 'rigid' and sum(self.h(s)) != sum(self.h(t)):
            return True
        return False

    def run(self):
        t0 = time.time()
        out = dict(name=self.name, d=self.d, w=self.w, n=self.n, S=self.S, nS=len(self.S), primes=self.primes,
                   moduli=[len(self.S) - self.K0[p] for p in self.primes], gmax=self.gmax,
                   perms=[''.join(map(str, pi)) for pi in self.perms], hodge_sym_ok=self.sym_ok, perm_ok=self.perm_ok)
        res = {}
        for T in ('S', 'C'):
            info, recs, conflicts = self.blocks(T)
            lb = sum(r['size'] * r['lb'] for r in recs)
            ub = sum(r['size'] * r['ub'] for r in recs)
            kinds = [v['kind'] for v in info.values()]
            cnt = {k: kinds.count(k) for k in ('rigid', 'pure', 'cm', 'res')}
            closed = all(r['closed'] for r in recs)
            noivhs = all(r['closed'] and not r['rule'].startswith('R6') and r['rule'] != 'R1:pure-res' for r in recs)
            rulecount, rulecount_f = {}, {}
            for r in recs:
                rulecount[r['rule']] = rulecount.get(r['rule'], 0) + r['size']
                rulecount_f[r['fixed'][0]] = rulecount_f.get(r['fixed'][0], 0) + r['size']
            fixed = dict(lb=sum(r['size'] * r['fixed'][1] for r in recs), ub=sum(r['size'] * r['fixed'][2] for r in recs),
                         closed=all(r['closed_fixed'] for r in recs),
                         closed_without_ivhs=all(r['closed_fixed'] and not r['fixed'][0].startswith('R6') and r['fixed'][0] != 'R1:pure-res' for r in recs),
                         rule_counts=rulecount_f)
            only = [dict(s=r['s'], t=r['t'], size=r['size'], hs=r['hs'], ht=r['ht'], zero_rules=r['zero_rules'],
                         zero_rules_fixed=r['zero_rules_fixed'], ub1_if_rule_dropped=r.get('ub1_if_rule_dropped'))
                    for r in recs if r.get('only_r5_or_grigid')]
            thmT = None
            if T == 'S' and cnt['res'] == 0:
                rho = {}
                for key in self.pieces:
                    if key[0] == 'S':
                        pass
                orbS = {}
                for O, v in info.items():
                    for x in O:
                        orbS[x] = v['kind']
                rho = {}
                for key in self.pieces:
                    if key[0] == 'S' and orbS[key] == 'pure':
                        rho[key[1]] = rho.get(key[1], 0) + sum(self.h(key))
                thmT = sum(v * v for v in rho.values()) + sum(len(self.aorbit(key)) for key in self.pieces if key[0] == 'S' and orbS[key] == 'rigid')
            res[T] = dict(lb=lb, ub=ub, closed=closed, closed_without_ivhs=noivhs, orbit_counts=cnt, rule_counts=rulecount,
                          fixed=fixed, only_r5_or_grigid=only,
                          conflicts=conflicts, thmT_formula=thmT,
                          orbits=[dict(v, key=[kstr(x) for x in O], cert=(kstr(v['cert']) if v.get('cert') else None)) for O, v in info.items()],
                          blocks=recs)
        out['Y'] = res['S']
        out['C'] = res['C']
        out['time'] = round(time.time() - t0, 1)
        return out


def fam_from_label(label, d):
    # cycN_abc or cycN_a-b-c
    n = int(label.split('_')[0][3:])
    ws = label.split('_')[1]
    w = tuple(int(x) for x in ws.split('-')) if '-' in ws else tuple(int(ch) for ch in ws)
    return n, w


def run_one(args):
    label, d = args
    n, w = fam_from_label(label, d)
    t0 = time.time()
    try:
        st = int(os.environ.get('V1_START', '50000'))
        F = Fam(label, d, w, n, start=st)
        out = F.run()
        fn = os.path.join(OUTD, f'v1_{label}_d{d}' + ('' if st == 50000 else f'_s{st}') + '.json')
        with open(fn, 'w') as f:
            json.dump(out, f, indent=1, default=str)
        Y, C = out['Y'], out['C']
        line = (f"{label:14s} d={d} n={n} w={w} |S|={out['nS']} mod={out['moduli']} gmax={out['gmax']} perms={out['perms']} "
                f"Y={Y['lb']}..{Y['ub']} closed={Y['closed']} noIVHS={Y['closed_without_ivhs']} Yorb={Y['orbit_counts']} thmT={Y['thmT_formula']} "
                f"C={C['lb']}..{C['ub']} closed={C['closed']} Corb={C['orbit_counts']} conflicts={len(Y['conflicts'])}+{len(C['conflicts'])} "
                f"FIXED: Y={Y['fixed']['lb']}..{Y['fixed']['ub']} closed={Y['fixed']['closed']} noIVHS={Y['fixed']['closed_without_ivhs']} "
                f"C={C['fixed']['lb']}..{C['fixed']['ub']} closed={C['fixed']['closed']} onlyR5/grigid(Y,C)={sum(o['size'] for o in Y['only_r5_or_grigid'])},{sum(o['size'] for o in C['only_r5_or_grigid'])} "
                f"sym={out['hodge_sym_ok']} perm_ok={out['perm_ok']} primes={out['primes']} t={out['time']}s")
    except Exception:
        line = f'{label} d={d} ERROR ' + traceback.format_exc().replace('\n', ' | ')
    return line


FAMS_A = [('cyc5_001', 6), ('cyc10_015', 6), ('cyc4_012', 7), ('cyc6_001', 7), ('cyc6_012', 7), ('cyc6_023', 7),
          ('cyc7_001', 7), ('cyc12_016', 7), ('cyc7_001', 8), ('cyc14_017', 8), ('cyc8_014', 9), ('cyc16_018', 9),
          ('cyc18_0-1-10', 9), ('cyc9_001', 10), ('cyc18_019', 10)]
FAMS_B = [(f'cyc{n}_001', d) for d in range(6, 11) for n in range(2, d + 1) if d % n == 0 or (d - 1) % n == 0]
FAMS_C = [('cyc8_012', 8)]
FAMS_D = [('cyc5_012', 6), ('cyc7_013', 7)]
FAMS_CTRL = [('cyc8_017', 8), ('cyc1_000', 6)]

if __name__ == '__main__':
    import multiprocessing as mp
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    fams = []
    for f in FAMS_A + FAMS_B + FAMS_C + FAMS_D + FAMS_CTRL:
        if f not in fams:
            fams.append(f)
    if which != 'all':
        fams = [f for f in fams if f'{f[0]}_d{f[1]}' in which.split(',')]
    # larger d / smaller n first (most expensive)
    fams.sort(key=lambda f: (-f[1] * 100 + fam_from_label(*f)[0]))
    summ = os.path.join(OUTD, f'v1_summary_{which if which == "all" else "sel"}' + os.environ.get('V1_TAG', '') + '.txt')
    with open(summ, 'w') as fo:
        fo.write(f'# V1 run started {time.ctime()} on {len(fams)} families\n')
    with mp.Pool(int(os.environ.get('V1_WORKERS', '8')), maxtasksperchild=1) as pool:
        for line in pool.imap_unordered(run_one, fams):
            with open(summ, 'a') as fo:
                fo.write(line + '\n')
            print(line, flush=True)
    with open(summ, 'a') as fo:
        fo.write(f'# finished {time.ctime()}\n')
