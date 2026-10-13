# plan: list (family, d, i) needed for order-2/3 tests of the open Y-blocks (Galois transport to a in {1,2});
# apply: read ord3 JSONs, ub := min(ub, ord3 'full' kernel dims over transported members, variants, primes).
# Only 'full' kernel dimensions are used (kernel-dim upper bounds, L6); projection ranks are NOT used.
import os, sys, json, glob
from math import gcd
sys.stdout.reconfigure(encoding='utf-8')
HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = os.path.join(os.path.dirname(HERE), 'outputs')
mode = sys.argv[1]


def units(N):
    return [k for k in range(1, N) if gcd(k, N) == 1]


def parse(s):
    a = s.split(':')
    return (a[0], int(a[1]), int(a[2]))


def lcm(a, b):
    return a * b // gcd(a, b)


allfams = sorted(glob.glob(os.path.join(OUTD, 'z1_cyc*_d*.json')))
tasks, summary = set(), []
for fn in allfams:
    D = json.load(open(fn))
    name, d, n, w = D['name'], D['d'], D['n'], D['w']
    sw = sum(w) % n
    N = lcm(d, n)
    ords = {}
    if mode == 'apply':
        for f2 in glob.glob(os.path.join(OUTD, f'ord3_{name}_d{d}_i*_p*_a8_o3.json')):
            R = json.load(open(f2))
            for b in R['blocks']:
                ords.setdefault((R['i'], b['sigma'], b['Delta']), []).append(
                    (R['p'], {lab: v['full'] for lab, v in b.items() if isinstance(v, dict) and 'full' in v}))
    nopen, nclosed, cov = 0, 0, 0
    lb_tot = ub_tot = 0
    newopen = []
    for r in D['blocks']:
        lb, ub = r['z'][1], r['z'][2]
        if r.get('open'):
            s, t = parse(r['s']), parse(r['t'])
            mem = set()
            for k in units(N):
                x, y = (k * s[1]) % d, (k * s[2]) % n
                yy = (k * t[2]) % n
                if x in (1, 2):
                    mem.add((x - 1, (y - sw) % n, (yy - y) % n))
            if mem:
                cov += 1
            for m in mem:
                tasks.add((name, d, m[0]))
            if mode == 'apply':
                best, where = None, None
                for m in sorted(mem):
                    for p, fulls in ords.get(m, []):
                        for lab, v in fulls.items():
                            if best is None or v < best:
                                best, where = v, (m, p, lab)
                r['ord3'] = dict(members=sorted(mem), best=best, where=where)
                if best is not None and best < ub:
                    ub = max(best, 0)
                    r['z_ord3'] = ['ord3:' + where[2], lb, ub]
                if lb == ub:
                    nclosed += r['size']
                else:
                    newopen.append(r)
            nopen += r['size']
        lb_tot += r['size'] * lb
        ub_tot += r['size'] * ub
        if lb > ub:
            print('  INCONSISTENT', name, d, r['s'], r['t'], lb, ub)
    if mode == 'apply':
        D['Ylb_ord3'], D['Yub_ord3'] = lb_tot, ub_tot
        D['open_after_ord3'] = newopen
        with open(fn.replace('z1_', 'z1o3_'), 'w') as f:
            json.dump(D, f, indent=1)
        print(f"{name:14s} d={d:2d} Y={D['Ylb']}..{D['Yub']} -> {lb_tot}..{ub_tot} closed={lb_tot == ub_tot} "
              f"open-before={nopen} closed-by-ord3={nclosed} block-orbits-coverable={cov}")
        for r in newopen:
            o = r.get('ord3', {})
            print(f"     open {r['s']}->{r['t']} size={r['size']} h={r['hs']}->{r['ht']} kinds={r['ks']},{r['kt']} "
                  f"aut={r['aut']} z={r.get('z_ord3', r['z'])} type={r['open']} ord3best={o.get('best')} at={o.get('where')}")
if mode == 'plan':
    print(' '.join(f'{a}:{b}:{c}' for a, b, c in sorted(tasks)))
