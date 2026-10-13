# Z1 triage of the residual symmetric families (Y-part only).
# Built on the Z1 copy of V1's code (v1lib.py, v1run.py in this folder). Changes against V1:
#  * Higgs-kernel rules disabled (r5valid, r5onesided -> False; V2 P1); galois_rigid then uses only Hodge numbers,
#    orbit sums and ub1 = 0.
#  * Irreducibility certificates PER MEMBER (embedding): V1 hull<=3 (+L4' for h20 = 2), hull to order 5 (+L4'),
#    and the full first-order commutant dim C(E,F) = 1 (U1 Prop. 5.1; mod-p kernel dim 1 => 1 in char 0).
#  * R2c uses every certified member (Galois transport of the pair), plus Cor. 5.2(iv),(v) per member.
#  * Eigen-split of pieces stabilised by a normaliser permutation pi, at a pi-symmetric point (where pi acts by
#    x -> x_pi): commutant algebra A, dims of e_b A e_a (Hom-commutants of eigen-subpieces), Hodge numbers of the
#    eigen-subpieces; blocks re-bounded (pure x pure all algebraic; irreducible vs pure 0; irreducible with
#    different Hodge numbers 0; otherwise Hom-commutant).
import os, sys, json, time, random, traceback
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v1lib import *
from v1run import Fam, kstr, fam_from_label, FAMS_Z

OUTD = os.path.join(os.path.dirname(HERE), 'outputs')


def comm_basis(thS, thT, hs, ht, p):
    """basis of {X = (X_l) : X_{l+1} thS_v = thT_v X_l for all v, l}; X_l stored row-major (ht[l] x hs[l])."""
    nl = len(hs)
    offs, U = [], 0
    for l in range(nl):
        offs.append(U)
        U += ht[l] * hs[l]
    if U == 0:
        return [], offs, 0
    nv = len(thS[0]) if nl > 1 else 0
    eqs = []
    for vi in range(nv):
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
        return [[1 if i == j else 0 for i in range(U)] for j in range(U)], offs, U
    X, nul = nm(eqs, U, p).nullspace()
    Xl = tl(X)
    return [[Xl[i][j] for i in range(U)] for j in range(nul)], offs, U


def matmul(A, B, p):
    return tl(nm(A, len(A[0]), p) * nm(B, len(B[0]), p))


def ispure_h(h):
    return h[0] == 0 and h[-1] == 0


class ZFam(Fam):
    HULLK = 5

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self._cdim, self._hull5, self._mcert, self._symR, self._split = {}, {}, {}, {}, {}
        self.Tk5 = None
        self.log = []

    # kernel rules off (V2 P1)
    def r5valid(self, s, t):
        return False

    def r5onesided(self, *a):
        return False

    def cdim(self, p, x):
        if (p, x) not in self._cdim:
            self._cdim[(p, x)] = ub1_block(self.theta(p, x), self.theta(p, x), self.h(x), self.h(x), p)
        return self._cdim[(p, x)]

    def hull5(self, p, x):
        if (p, x) not in self._hull5:
            if self.Tk5 is None:
                Tk = [set(t) for t in self.Tk]
                while len(Tk) <= self.HULLK:
                    Tk.append({madd(a, b) for a in Tk[-1] for b in self.S})
                self.Tk5 = Tk
            P = self.pieces[x]
            self._hull5[(p, x)] = self.rings[p].hull(P['shift'], P['c'], self.S, self.Tk5, P['qmax'], kmax=self.HULLK)[0]
        return self._hull5[(p, x)]

    def member_certs(self, x):
        if x in self._mcert:
            return self._mcert[x]
        h = self.h(x)
        out = []
        if h[0] in (1, 2):
            for p in self.primes:
                hr = self.hull(p, x)
                if hr[-1] == sum(h) and (h[0] == 1 or self.l4(p, x)['cert']):
                    out.append('hull3' if h[0] == 1 else 'hull3+L4p')
                    break
        for p in self.primes:
            if self.cdim(p, x) == 1:
                out.append('comm')
                break
        if not out and h[0] in (1, 2):
            for p in self.primes:
                hr = self.hull5(p, x)
                if hr[-1] == sum(h) and (h[0] == 1 or self.l4(p, x)['cert']):
                    out.append('hull5' if h[0] == 1 else 'hull5+L4p')
                    break
        self._mcert[x] = out
        return out

    def classify_orbits(self, T):
        pure = 'pure' if T == 'S' else 'cm'
        keys = sorted(k for k in self.pieces if k[0] == T)
        orbits, self_orb = [], {}
        for key in keys:
            if key in self_orb:
                continue
            O = self.gorbit(key)
            for x in O:
                self_orb[x] = O
            orbits.append(O)
        info = {}
        for O in orbits:
            ispure = all(ispure_h(self.h(x)) for x in O) if T == 'S' else all(self.h(x)[0] * self.h(x)[1] == 0 for x in O)
            rec = dict(members=[kstr(x) for x in O], h={kstr(x): self.h(x) for x in O})
            if ispure:
                rec['kind'] = pure
            else:
                certs = {x: self.member_certs(x) for x in O}
                rec['member_certs'] = {kstr(x): certs[x] for x in O if certs[x]}
                rec['cdims'] = {kstr(x): [self.cdim(p, x) for p in self.primes] for x in O}
                if self.h(O[0])[0] in (1, 2) or any(self.h(x)[0] in (1, 2) for x in O):
                    rec['hull3'] = {kstr(x): self.hull(self.primes[0], x) for x in O if self.h(x)[0] in (1, 2)}
                    h5 = {kstr(x): self.hull5(self.primes[0], x) for x in O if (x, self.primes[0]) and (self.primes[0], x) in self._hull5}
                    if h5:
                        rec['hull5'] = h5
                cm = [x for x in O if certs[x]]
                rec['certs_all'] = cm
                rec['kind'] = 'rigid' if cm else 'res'
                rec['cert'] = cm[0] if cm else None
            info[O] = rec
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
        self._info = info
        return orbits, self_orb, info

    def rule2c(self, s, t, Ks, Kt, info, orb):
        T = s[0]
        nl = len(self.h(s))
        cs = set(info[orb[s]].get('certs_all') or [])
        ct = set(info[orb[t]].get('certs_all') or [])
        for k in self.units(T):
            s2, t2 = self.act(k, s), self.act(k, t)
            hs2, ht2 = self.h(s2), self.h(t2)
            if s2 in cs and any(hs2[l] > ht2[l] for l in range(nl)):
                return True
            if t2 in ct and any(hs2[l] < ht2[l] for l in range(nl)):
                return True
            if s2 in cs and t2 in ct and hs2 != ht2:
                return True
            if T == 'S' and s2 in cs and ispure_h(ht2) and not ispure_h(hs2):
                return True
            if T == 'S' and t2 in ct and ispure_h(hs2) and not ispure_h(ht2):
                return True
        if Ks == 'rigid' and Kt == 'rigid' and sum(self.h(s)) != sum(self.h(t)):
            return True
        return False

    # ---------------- eigen-split at a pi-symmetric point
    def stab(self, x):
        return [pi for pi, mp in zip(self.perms, self.permmaps) if mp.get(x) == x]

    def symring(self, p):
        if p not in self._symR:
            rep = {}
            for m in self.S:
                if m in rep:
                    continue
                O, todo = {m}, [m]
                while todo:
                    z = todo.pop()
                    for pi in self.perms:
                        y = perm_apply(pi, z)
                        if y not in O:
                            O.add(y)
                            todo.append(y)
                for z in O:
                    rep[z] = min(O)
            rng = random.Random(f'Z1-sym|{self.name}|{self.d}|{p}')
            val, F = {}, {}
            for m in self.S:
                if rep[m] not in val:
                    val[rep[m]] = rng.randrange(1, p)
                F[m] = val[rep[m]]
            R = Ring(self.d, self.w, self.n, p, F)
            bad = R.check_hilbert()
            if bad:
                raise RuntimeError(f'symmetric point singular at p={p}')
            self._symR[p] = R
        return self._symR[p]

    def split_data(self, p, x):
        if (p, x) in self._split:
            return self._split[(p, x)]
        R = self.symring(p)
        P = self.pieces[x]
        d, h = self.d, self.h(x)
        th = R.theta(P['shift'], P['c'], self.S, P['qmax'])
        st = self.stab(x)
        rec = dict(th=th, parts=[('all', None, tuple(h))], pi=None)
        if st:
            pi = st[0]
            Pis = []
            for q in range(1, P['qmax'] + 1):
                std = R.std(q * d - P['shift'], P['c'])
                cols = []
                for M in std:
                    v, dims = R.reduce_poly(P['shift'], q, P['c'], {perm_apply(pi, M): 1})
                    off = sum(dims[:-1])
                    cols.append(v[off:off + dims[-1]])
                hq = len(std)
                Pis.append([[cols[j][i] for j in range(hq)] for i in range(hq)])
            inv2 = pow(2, p - 2, p)
            invol = all(matmul(M, M, p) == [[1 if i == j else 0 for j in range(len(M))] for i in range(len(M))] for M in Pis if M)
            comm_ok = True
            for l in range(len(Pis) - 1):
                if not Pis[l] or not Pis[l + 1]:
                    continue
                for M in th[l]:
                    if matmul(Pis[l + 1], M, p) != matmul(M, Pis[l], p):
                        comm_ok = False
            rec['pi'] = ''.join(map(str, pi))
            rec['invol'], rec['pi_commutes_theta'] = invol, comm_ok
            if invol and comm_ok:
                parts = []
                for sg in (1, -1):
                    E = [[[((1 if i == j else 0) + sg * M[i][j]) * inv2 % p for j in range(len(M))] for i in range(len(M))] if M else []
                         for M in Pis]
                    hh = tuple(rank_rows(Em, len(Em), p) if Em else 0 for Em in E)
                    if sum(hh):
                        parts.append(('+' if sg == 1 else '-', E, hh))
                rec['parts'] = parts
        self._split[(p, x)] = rec
        return rec

    def split_pair(self, p, x, y):
        Dx, Dy = self.split_data(p, x), self.split_data(p, y)
        hx, hy = self.h(x), self.h(y)
        basis, offs, U = comm_basis(Dx['th'], Dy['th'], hx, hy, p)
        out = {}
        for (la, Ea, ha) in Dx['parts']:
            for (lb, Eb, hb) in Dy['parts']:
                vecs = []
                for X in basis:
                    vec = []
                    for l in range(len(hx)):
                        if hx[l] == 0 or hy[l] == 0:
                            continue
                        Xl = [X[offs[l] + i * hx[l]: offs[l] + (i + 1) * hx[l]] for i in range(hy[l])]
                        if Ea is not None:
                            Xl = matmul(Xl, Ea[l], p)
                        if Eb is not None:
                            Xl = matmul(Eb[l], Xl, p)
                        vec += [z for r in Xl for z in r]
                    vecs.append(vec)
                out[(la, lb)] = rank_rows(vecs, len(vecs[0]), p) if vecs and vecs[0] else 0
        return out, Dx, Dy

    def part_pure_orbit(self, p, x, la):
        for k in self.units('S'):
            y = self.act(k, x)
            if la == 'all':
                hh = [self.h(y)]
            else:
                hh = [h for (l, E, h) in self.split_data(p, y)['parts'] if l == la]
            if hh and not ispure_h(hh[0]):
                return False
        return True

    def split_block(self, B, certs_all):
        """bounds (lb, ub) for the block orbit B from the eigen-split; min ub / max lb over members and primes."""
        best_ub, best_lb, detail = None, 0, []
        for (x, y) in B:
            for p in self.primes:
                dims, Dx, Dy = self.split_pair(p, x, y)
                selfx = self.split_pair(p, x, x)[0] if x != y else dims
                selfy = self.split_pair(p, y, y)[0] if x != y else dims
                irr = lambda D, la, sd, z: (la == 'all' and z in certs_all) or (sd.get((la, la)) == 1)
                lb = ub = 0
                parts_desc = []
                for (la, Ea, ha) in Dx['parts']:
                    for (lbb, Eb, hb) in Dy['parts']:
                        pa, pb = ispure_h(ha), ispure_h(hb)
                        ia, ib = irr(Dx, la, selfx, x), irr(Dy, lbb, selfy, y)
                        if pa and pb:
                            v = sum(ha) * sum(hb)
                            ub += v
                            # algebraic (NS (x) NS) only if both parts are pure at EVERY Galois embedding
                            if self.part_pure_orbit(p, x, la) and self.part_pure_orbit(p, y, lbb):
                                lb += v
                                why = 'pure-pure(orbit):alg'
                            else:
                                why = 'pure-pure(embedding only):ub'
                        elif (ia and not pa and pb) or (ib and not pb and pa):
                            why = 'irr-pure:0'
                        elif ia and ib and ha != hb:
                            why = 'irr-hodge:0'
                        else:
                            v = dims[(la, lbb)]
                            ub += v
                            if x == y and la == lbb:
                                lb += 1
                            why = f'comm={v}'
                        parts_desc.append((la, ha, lbb, hb, why))
                if best_ub is None or ub < best_ub:
                    best_ub = ub
                best_lb = max(best_lb, lb)
                detail.append(dict(x=kstr(x), y=kstr(y), p=p, lb=lb, ub=ub, pi_x=Dx.get('pi'), pi_y=Dy.get('pi'),
                                   checks=[Dx.get('invol'), Dx.get('pi_commutes_theta'), Dy.get('invol'), Dy.get('pi_commutes_theta')],
                                   parts=parts_desc))
        return best_lb, best_ub, detail

    def run_z(self):
        t0 = time.time()
        info, recs, conflicts = self.blocks('S')
        orb = {}
        for O in info:
            for x in O:
                orb[x] = O
        certs_all = set()
        for O, v in info.items():
            certs_all |= set(v.get('certs_all') or [])
        parse = lambda s: (s.split(':')[0], int(s.split(':')[1]), int(s.split(':')[2]))
        for r in recs:
            rule, lb, ub = r['fixed']
            r['z'] = [rule, lb, ub]
            if lb == ub or not self.perms:
                continue
            s, t = parse(r['s']), parse(r['t'])
            B = sorted(set((self.act(k, s), self.act(k, t)) for k in self.units('S')))
            if not any(self.stab(x) or self.stab(y) for x, y in B):
                continue
            try:
                slb, sub, det = self.split_block(B, certs_all)
                r['split'] = dict(lb=slb, ub=sub, detail=det)
                nlb, nub = max(lb, slb), min(ub, sub)
                r['z'] = ['split' if nlb == nub else rule + '+split', nlb, nub]
            except Exception:
                r['split_error'] = traceback.format_exc()
        for r in recs:
            ks, kt = info[orb[parse(r['s'])]]['kind'], info[orb[parse(r['t'])]]['kind']
            if r['z'][1] == r['z'][2]:
                r['open'] = None
            elif 'res' in (ks, kt):
                r['open'] = 'aut-nonrigid' if r['aut'] else 'undecided-irreducibility'
            elif r['aut']:
                r['open'] = 'aut-rigid'
            else:
                r['open'] = 'equal-invariants' if r['hs'] == r['ht'] else 'other'
        lb = sum(r['size'] * r['z'][1] for r in recs)
        ub = sum(r['size'] * r['z'][2] for r in recs)
        kinds = [v['kind'] for v in info.values()]
        certcount = {}
        for v in info.values():
            for x, cl in (v.get('member_certs') or {}).items():
                for c in cl:
                    certcount[c] = certcount.get(c, 0) + 1
        out = dict(name=self.name, d=self.d, w=self.w, n=self.n, S=self.S, nS=len(self.S), primes=self.primes,
                   perms=[''.join(map(str, pi)) for pi in self.perms], Ylb=lb, Yub=ub, closed=(lb == ub),
                   orbit_counts={k: kinds.count(k) for k in ('rigid', 'pure', 'res')}, cert_counts=certcount,
                   conflicts=len(conflicts),
                   rule_counts={}, open_blocks=[r for r in recs if r['open']],
                   orbits=[dict({k: v for k, v in rec.items() if k not in ('cert', 'certs_all')}, key=[kstr(x) for x in O],
                                cert=kstr(rec['cert']) if rec.get('cert') else None,
                                certs_all=[kstr(x) for x in rec.get('certs_all') or []]) for O, rec in info.items()],
                   blocks=recs, time=round(time.time() - t0, 1))
        for r in recs:
            out['rule_counts'][r['z'][0]] = out['rule_counts'].get(r['z'][0], 0) + r['size']
        out['inconsistent'] = sum(r['size'] for r in recs if r['z'][1] > r['z'][2])
        out['conflicts'] = f"{len(conflicts)} incons={out['inconsistent']}"
        return out


def run_one(args):
    label, d = args
    n, w = fam_from_label(label, d)
    try:
        F = ZFam(label, d, w, n, start=50000)
        out = F.run_z()
        with open(os.path.join(OUTD, f'z1_{label}_d{d}.json'), 'w') as f:
            json.dump(out, f, indent=1, default=str)
        ob = out['open_blocks']
        reasons = {}
        for r in ob:
            reasons[r['open']] = reasons.get(r['open'], 0) + r['size']
        line = (f"{label:14s} d={d} |S|={out['nS']} perms={out['perms']} Y={out['Ylb']}..{out['Yub']} closed={out['closed']} "
                f"orb={out['orbit_counts']} certs={out['cert_counts']} open={reasons} conflicts={out['conflicts']} "
                f"rules={out['rule_counts']} t={out['time']}s")
    except Exception:
        line = f'{label} d={d} ERROR ' + traceback.format_exc().replace('\n', ' | ')
    return line


if __name__ == '__main__':
    import multiprocessing as mp
    sys.stdout.reconfigure(encoding='utf-8')
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    fams = list(FAMS_Z)
    if which != 'all':
        fams = [f for f in fams if f'{f[0]}_d{f[1]}' in which.split(',')]
    fams.sort(key=lambda f: -f[1])
    summ = os.path.join(OUTD, 'z1_summary' + ('' if which == 'all' else '_sel') + '.txt')
    with open(summ, 'w') as fo:
        fo.write(f'# Z1 triage started {time.ctime()} on {len(fams)} families\n')
    with mp.Pool(int(os.environ.get('Z1_WORKERS', '8')), maxtasksperchild=1) as pool:
        for line in pool.imap_unordered(run_one, fams):
            with open(summ, 'a') as fo:
                fo.write(line + '\n')
            print(line, flush=True)
    with open(summ, 'a') as fo:
        fo.write(f'# finished {time.ctime()}\n')
