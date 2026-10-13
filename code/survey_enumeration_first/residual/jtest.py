# Z1 J-test: higher-order (k <= KMAX) F^2-jet test for open Y-blocks, any h^{0,2} (complements ord3, which needs h02 = 0).
# A flat (0,0) section phi of Hom(E_s, E_t) satisfies phi(F^2 E_s(t)) in F^2 E_t(t) for all t. Along an arc F + tL:
#   omega_P(t) = sum_k W_k(P) t^k,  W_k(P) = (-1)^k eta(P L^k, k+1)  (flat Taylor coefficients in the fibre at F),
#   phi W^s_k(P) = sum_{m<=k} W^t_{k-m}(c^(m) P),  c^(0) shared by all arcs, c^(m) (m >= 1) per arc.
# The kernel dimension of this linear system (unknowns phi, c) is >= dim of flat (0,0) sections (phi |-> (phi, c) is
# injective), and a mod-p kernel dimension is an upper bound for char 0 (L6). Only kernel dimensions are reported.
import os, sys, json, time, random, traceback
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v1lib import *
from v1run import Fam, kstr, fam_from_label

OUTD = os.path.join(os.path.dirname(HERE), 'outputs')
KMAX, NARCS, MAXMEM = 4, 8, 4


def ppow_mul(P, L, p):
    out = {}
    for e1, c1 in P.items():
        for e2, c2 in L.items():
            e = madd(e1, e2)
            out[e] = (out.get(e, 0) + c1 * c2) % p
    return {e: c for e, c in out.items() if c}


def wcoef(fam, p, x, L, kmax):
    """W[k][P] (vectors of length dim E_x) for k = 0..kmax, P over the level-1 basis of E_x."""
    R = fam.rings[p]
    Pc = fam.pieces[x]
    sh, c = Pc['shift'], Pc['c']
    dimE = sum(fam.h(x))
    top = R.cols(fam.d - sh, c)[0]
    W = []
    polys = [{P: 1} for P in top]
    for k in range(kmax + 1):
        row = []
        for poly in polys:
            v, dims = R.reduce_poly(sh, k + 1, c, poly)
            v = (v + [0] * dimE)[:dimE] if len(v) < dimE else v[:dimE]
            sg = 1 if k % 2 == 0 else p - 1
            row.append([(sg * z) % p for z in v])
        W.append(row)
        polys = [ppow_mul(poly, L, p) for poly in polys]
    return W


def jdims(fam, p, x, y, arcs):
    ds, dt = sum(fam.h(x)), sum(fam.h(y))
    rs, rt = fam.h(x)[0], fam.h(y)[0]
    Ws = [wcoef(fam, p, x, L, KMAX) for L in arcs]
    Wt = [wcoef(fam, p, y, L, KMAX) for L in arcs]
    nphi, nc = ds * dt, rs * rt
    N = nphi + nc * (1 + len(arcs) * KMAX)
    iphi = lambda i, j: i * ds + j
    ic = lambda a, m, q, pp: nphi + (q * rs + pp if m == 0 else nc * (1 + a * KMAX + (m - 1)) + q * rs + pp)
    rows, out = [], []
    for k in range(KMAX + 1):
        for a in range(len(arcs)):
            for pp in range(rs):
                for i in range(dt):
                    row = [0] * N
                    for j in range(ds):
                        v = Ws[a][k][pp][j]
                        if v:
                            row[iphi(i, j)] = v
                    for m in range(k + 1):
                        for q in range(rt):
                            v = Wt[a][k - m][q][i]
                            if v:
                                idx = ic(a, m, q, pp)
                                row[idx] = (row[idx] - v) % p
                    if any(row):
                        rows.append(row)
        # unknowns c^(m) with m > k are not yet constrained: count them out
        free_later = nc * len(arcs) * (KMAX - k)
        rk = rank_rows(rows, N, p)
        out.append(N - rk - free_later)
    return out   # out[k] = kernel dim using orders 0..k


def run_family(tag):
    fn = os.path.join(OUTD, f'z1_{tag}.json')
    D = json.load(open(fn))
    label, d = D['name'], D['d']
    n, w = fam_from_label(label, d)
    t0 = time.time()
    F = Fam(label, d, w, n, start=50000)
    res = []
    parse = lambda s: (s.split(':')[0], int(s.split(':')[1]), int(s.split(':')[2]))
    for r in D['blocks']:
        if not r.get('open'):
            continue
        s, t = parse(r['s']), parse(r['t'])
        B = sorted(set((F.act(k, s), F.act(k, t)) for k in F.units('S')))
        mem = [(x, y) for x, y in B if F.h(x)[0] >= 1 and F.h(x)[0] == F.h(y)[0]]
        mem.sort(key=lambda xy: sum(F.h(xy[0])) * sum(F.h(xy[1])))
        rec = dict(s=r['s'], t=r['t'], size=r['size'], hs=r['hs'], ht=r['ht'], open=r['open'], z=r['z'], tests=[])
        best = None
        for (x, y) in mem[:MAXMEM]:
            for p in F.primes:
                rng = random.Random(f'Z1-J|{label}|{d}|{p}|{kstr(x)}|{kstr(y)}')
                arcs = [{m: rng.randrange(p) for m in F.S} for _ in range(NARCS)]
                dims = jdims(F, p, x, y, arcs)
                rec['tests'].append(dict(x=kstr(x), y=kstr(y), hx=F.h(x), hy=F.h(y), p=p, dims=dims))
                best = min(dims) if best is None else min(best, min(dims))
        rec['jbest'] = best
        res.append(rec)
    out = dict(name=label, d=d, KMAX=KMAX, NARCS=NARCS, blocks=res, time=round(time.time() - t0, 1))
    with open(os.path.join(OUTD, f'jt_{tag}.json'), 'w') as f:
        json.dump(out, f, indent=1)
    lines = [f'{tag} t={out["time"]}s']
    for rec in res:
        lines.append(f"   {rec['s']}->{rec['t']} size={rec['size']} h={rec['hs']} open={rec['open']} z={rec['z']} "
                     f"jbest={rec['jbest']} dims={[(tt['x'], tt['y'], tt['p'] % 1000, tt['dims']) for tt in rec['tests'][:2]]}")
    return '\n'.join(lines)


def safe(tag):
    try:
        return run_family(tag)
    except Exception:
        return f'{tag} ERROR ' + traceback.format_exc().replace('\n', ' | ')


if __name__ == '__main__':
    import multiprocessing as mp
    sys.stdout.reconfigure(encoding='utf-8')
    tags = sys.argv[1:]
    with mp.Pool(min(8, len(tags)), maxtasksperchild=1) as pool:
        for line in pool.imap_unordered(safe, tags):
            print(line, flush=True)
