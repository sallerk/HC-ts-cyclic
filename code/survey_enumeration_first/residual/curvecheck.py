# Curve explanation of hull deficiencies. If F|_{x_i=0} is a single monomial M = x_j^al x_k^be (degree d), then
# Y_F cap {x_i = 0} = {M + u^d = 0} splits into g = gcd(d, al, be) curves permuted by mu_d (u-action), G-invariant;
# their classes give (at most) one class in each piece (a, chi=0) with a = 0 mod d/g, a != 0.
# Prediction: delta_pred(a) = #{i : F|_{x_i=0} monomial with a = 0 mod d/g_i}. Compare with the computed hull
# deficiency dim E - hullrank (order <= 5, lower bound for the true hull rank) for members with h20 >= 1.
import os, sys, json, glob
from math import gcd
sys.stdout.reconfigure(encoding='utf-8')
OUTD = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'outputs')
out = []
for fn in sorted(glob.glob(os.path.join(OUTD, 'z1_cyc*_d*.json'))):
    D = json.load(open(fn))
    d, S = D['d'], [tuple(m) for m in D['S']]
    tag = f"{D['name']}_d{d}"
    planes = []
    for i in range(3):
        R = [m for m in S if m[i] == 0]
        if len(R) == 1:
            m = R[0]
            al, be = [m[j] for j in range(3) if j != i]
            g = gcd(d, gcd(al, be))
            if g > 1:
                planes.append((i, m, g))
    pred = lambda a: sum(1 for (i, m, g) in planes if a % (d // g) == 0)
    rows, bad = [], 0
    for o in D['orbits']:
        for k in ('hull5', 'hull3'):
            for x, hr in (o.get(k) or {}).items():
                a, chi = int(x.split(':')[1]), int(x.split(':')[2])
                h = o['h'][x]
                deff = sum(h) - hr[-1]
                pd = pred(a) if chi == 0 else 0
                ok = deff == pd
                bad += not ok
                if deff or pd:
                    rows.append(f"{x} h={h} hull={hr[-1]} deficiency={deff} predicted={pd} {'OK' if ok else 'MISMATCH'} kind={o['kind']}")
            if o.get('hull5') and k == 'hull5':
                break
    out.append(f"{tag:18s} split planes (i, monomial, #curves)={planes} mismatches={bad}")
    seen = set()
    for r in rows:
        if r not in seen:
            out.append('    ' + r)
            seen.add(r)
open(os.path.join(OUTD, 'curvecheck.txt'), 'w', encoding='utf-8').write('\n'.join(out) + '\n')
print('\n'.join(out))
