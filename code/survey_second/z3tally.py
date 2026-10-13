# z3tally.py -- tabulate Z3 results (own outputs only). python -I -B z3tally.py
import sys, os, json, glob, collections
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z3run import FAMILIES, CONTROLS

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'outputs')

def load(name, d, suf=''):
    fn = os.path.join(OUT, f'z3_{name}_d{d}{suf}.json')
    if not os.path.exists(fn):
        return None
    with open(fn, encoding='utf-8') as f:
        return json.load(f)

lines = []
tot_tools = collections.Counter()
lines.append('family          d  n  w        |S| primes          Y(both)     Y(p0 only)  Y(p1 only)  #borb  certE certU planes perms(valid)  tools closing (block orbits)')
for (name, d) in FAMILIES + CONTROLS:
    r = load(name, d)
    r0, r1 = load(name, d, '_only_p0'), load(name, d, '_only_p1')
    if r is None:
        lines.append(f'{name} d{d}: MISSING'); continue
    tc = collections.Counter()
    for b in r['blocks']:
        t = (b['tool'] or 'none').split('@')[0]
        tc[t] += 1
        tot_tools[t] += 1
    nE = sum(1 for v in r['pieces'].values() if v['certE'] and not v['pure'])
    nU = sum(1 for v in r['pieces'].values() if v['certU'])
    pl = ';'.join(f"i{p['i']}:M{tuple(p['M'])}g{p['g']}" for p in r['planes']) or '-'
    pv = ';'.join(f"{tuple(p['pi'])}{'V' if p['ii'] and p['iii'] else 'x'}" for p in r['perms']) or '-'
    def ys(q):
        return '-' if q is None else (f"{q['Ylb']}" if q['Ylb'] == q['Yub'] else f"{q['Ylb']}..{q['Yub']}")
    lines.append(f"{name:14s} {d:2d} {r['n']:2d} {str(tuple(r['w'])):9s} {r['S']:3d} {str(r['primes']):15s} {ys(r):11s} {ys(r0):11s} {ys(r1):11s} "
                 f"{len(r['blocks']):5d} {nE:5d} {nU:5d} {pl} {pv}  " + ', '.join(f'{k}:{v}' for k, v in sorted(tc.items())))
    if r['errors']:
        lines.append(f'   ERRORS: {r["errors"]}')
    if r0 and r0['errors'] or r1 and r1['errors']:
        lines.append(f'   ERRORS single-prime: {r0 and r0["errors"]} {r1 and r1["errors"]}')
    for b in r['blocks']:
        if not b['closed']:
            lines.append(f"   OPEN block {tuple(b['rep'][0])}->{tuple(b['rep'][1])} x{b['size']} h={b['h']} lb={b['lb']} ub={b['ub']} {b['tool']} :: {' '.join(b['details'])[:600]}")
lines.append('')
lines.append('tool totals over all families (block orbits): ' + ', '.join(f'{k}:{v}' for k, v in sorted(tot_tools.items())))
# blocks closed by non-structural tools (homcomm, J, J_V, split): list
lines.append('')
lines.append('Block orbits closed by a computed bound (homcomm / J / J_V / split), per family:')
for (name, d) in FAMILIES + CONTROLS:
    r = load(name, d)
    if r is None:
        continue
    for b in r['blocks']:
        t = (b['tool'] or '')
        if t.startswith(('homcomm', 'J', 'split')):
            lines.append(f"   {name} d{d} {tuple(b['rep'][0])}->{tuple(b['rep'][1])} x{b['size']} lb={b['lb']} ub={b['ub']} {t}")
txt = '\n'.join(lines)
print(txt)
with open(os.path.join(OUT, 'z3_table.txt'), 'w', encoding='utf-8') as f:
    f.write(txt + '\n')
