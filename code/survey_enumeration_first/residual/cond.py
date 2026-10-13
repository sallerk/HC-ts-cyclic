# (1) per family: orbits whose irreducibility certificate needs the commutant (no hull3/hull5 member) or hull5;
# (2) for open blocks with a hull-deficient SOURCE x (h02(x) = 0): J_V = J_raw - (dim E_x - r_x) * dim E_y, where
#     r_x = computed hull rank (order 5, a LOWER bound for the true hull rank). J_V bounds the block only under the
#     inference "true hull = computed span" (then E_x = U + U^perp, U^perp pure flat; y irreducible non-pure kills
#     Hom(U^perp, E_y)).  Reported as CONDITIONAL.
import os, sys, json, glob
sys.stdout.reconfigure(encoding='utf-8')
OUTD = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'outputs')
out = []
for fn in sorted(glob.glob(os.path.join(OUTD, 'z1o3_cyc*_d*.json'))):
    D = json.load(open(fn))
    tag = f"{D['name']}_d{D['d']}"
    commonly, h5only, hull, certs = 0, 0, {}, set()
    for o in D['orbits']:
        if o['kind'] != 'rigid':
            pass
        mc = o.get('member_certs') or {}
        allc = ' '.join(' '.join(v) for v in mc.values())
        if o['kind'] == 'rigid':
            if 'hull3' not in allc and 'hull5' not in allc:
                commonly += 1
            elif 'hull3' not in allc and 'comm' not in allc:
                h5only += 1
        certs |= set(o.get('certs_all') or [])
        for k in ('hull5', 'hull3'):
            for x, v in (o.get(k) or {}).items():
                hull.setdefault(x, v[-1])
    out.append(f"{tag:18s} rigid orbits certified only by commutant: {commonly}; only by hull5: {h5only}")
    jf = os.path.join(OUTD, f'jt_{tag}.json')
    if not os.path.exists(jf):
        continue
    J = {(b['s'], b['t']): b for b in json.load(open(jf))['blocks']}
    for r in D.get('open_after_ord3', []):
        b = J.get((r['s'], r['t']))
        if not b or r['ks'] != 'res' or r['s'] == r['t']:
            continue
        best = None
        for t in b['tests']:
            x, y, hx, hy = t['x'], t['y'], t['hx'], t['hy']
            if hx[2] != 0 or x not in hull:
                continue
            delta = sum(hx) - hull[x]
            jv = min(t['dims']) - delta * sum(hy)
            ycert = y in certs
            rec = (jv, x, y, hx, hy, hull[x], ycert)
            if ycert and (best is None or jv < best[0]):
                best = rec
        if best:
            out.append(f"   COND {r['s']}->{r['t']} size={r['size']} type={r['open']}: J_V={best[0]} at {best[1]}->{best[2]} "
                       f"h={best[3]}->{best[4]} hullrank={best[5]} target certified={best[6]}")
open(os.path.join(OUTD, 'cond.txt'), 'w', encoding='utf-8').write('\n'.join(out) + '\n')
print('\n'.join(out))
