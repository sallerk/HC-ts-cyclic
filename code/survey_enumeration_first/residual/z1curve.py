# Curve-split closure of the residual chi=0 pieces (see curvecheck.py): E_x = U + P, P = span of the a-components of
# the classes of the curves in Y_F cap {x_i = 0} (algebraic, flat, pure, rank delta), U = P^perp = flat hull of F^2
# (rank U <= dim E - delta from P, >= dim E - delta from the computed hull: equality). U irreducible by L2 at members
# with h20 = 1 (Galois transport of block dimensions otherwise). Block rules:
#   (x,x): 1 + delta^2 (1_U and End(P) = NS (x) NS);  (x,t), (t,x) with t in a pure orbit: delta * dim E_t (NS (x) NS);
#   (x,y), y certified irreducible non-pure: Hom(P,E_y) = 0 and Hom(U,E_y) <= J_V = J_raw - delta * dim E_y;
#   (y,x), y certified irreducible: Hom(E_y,P) = 0, Hom(E_y,U) = 0 if h(E_y) != h(U) = h(E_x) - (0,delta,0).
import os, sys, json, glob
from math import gcd
sys.stdout.reconfigure(encoding='utf-8')
OUTD = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'outputs')
out = []
for fn in sorted(glob.glob(os.path.join(OUTD, 'z1o3_cyc*_d*.json'))):
    D = json.load(open(fn))
    d, S = D['d'], [tuple(m) for m in D['S']]
    tag = f"{D['name']}_d{d}"
    planes = []
    for i in range(3):
        R = [m for m in S if m[i] == 0]
        if len(R) == 1:
            al, be = [R[0][j] for j in range(3) if j != i]
            g = gcd(d, gcd(al, be))
            if g > 1:
                planes.append(g)
    delta = lambda a: sum(1 for g in planes if a % (d // g) == 0)
    kind, certs, h1member = {}, set(), {}
    for o in D['orbits']:
        for x in o['members']:
            kind[x] = o['kind']
        certs |= set(o.get('certs_all') or [])
        if o['kind'] == 'res':
            h1member[tuple(o['members'])] = any(o['h'][x][0] == 1 and (x in (o.get('hull5') or {}) or x in (o.get('hull3') or {})) for x in o['members'])
    resorb_ok = {x: v for k, v in h1member.items() for x in k}
    jf = os.path.join(OUTD, f'jt_{tag}.json')
    J = {(b['s'], b['t']): b for b in json.load(open(jf))['blocks']} if os.path.exists(jf) else {}
    isR = lambda x: kind[x] == 'res' and int(x.split(':')[2]) == 0 and delta(int(x.split(':')[1])) > 0
    Y = 0
    left, changed = [], []
    for r in D['blocks']:
        lb, ub = r['z'][1], r['z'][2]
        if 'z_ord3' in r:
            ub = min(ub, r['z_ord3'][2])
        jb = J.get((r['s'], r['t']))
        if jb and jb.get('jbest') is not None:
            ub = min(ub, jb['jbest'])
        if lb < ub:
            s, t = r['s'], r['t']
            a = int(s.split(':')[1])
            new = None
            if isR(s) and s == t and resorb_ok.get(s):
                new = 1 + delta(a) ** 2
            elif isR(s) and kind[t] == 'pure':
                new = delta(a) * sum(r['ht'])
            elif isR(t) and kind[s] == 'pure':
                new = delta(a) * sum(r['hs'])
            elif isR(s) and kind[t] == 'rigid' and jb:
                jv = [min(tt['dims']) - delta(a) * sum(tt['hy']) for tt in jb['tests'] if tt['y'] in certs and tt['hx'][2] == 0]
                jv += [min(tt['dims']) - delta(a) * sum(tt['hy']) for tt in jb['tests'] if tt['y'] in certs and tt['hx'][2] > 0]
                if jv and min(jv) == 0:
                    new = 0
            elif isR(t) and kind[s] == 'rigid':
                hU = [r['ht'][0], r['ht'][1] - delta(a), r['ht'][2]]
                if s in certs and list(r['hs']) != hU:
                    new = 0
            if new is not None and lb <= new <= ub:
                changed.append(f"    {s}->{t} size={r['size']} [{lb},{ub}] -> {new}")
                lb = ub = new
            elif new is not None:
                left.append(f"    INCONSISTENT {s}->{t} [{lb},{ub}] new={new}")
            else:
                left.append(f"    open {s}->{t} size={r['size']} h={r['hs']}->{r['ht']} kinds={kind[s]},{kind[t]} [{lb},{ub}] type={r.get('open')}")
        Y += r['size'] * lb
    status = 'CLOSED' if not left else 'open'
    out.append(f"{tag:18s} planes={planes} Y_lb={Y} {status} (changed {len(changed)} block orbits)")
    out += changed + left
open(os.path.join(OUTD, 'z1curve.txt'), 'w', encoding='utf-8').write('\n'.join(out) + '\n')
print('\n'.join(out))
