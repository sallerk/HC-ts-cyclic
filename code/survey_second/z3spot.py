# z3spot.py -- Z3 values at the spot checks quoted in TOOLS.txt sec. S, and route differences vs book.txt.
import sys, os, re, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf-8')
from z3lib import *
from z3run import PrimeData

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'outputs')
lines = []
def say(s):
    lines.append(s); print(s)
def load(tag):
    with open(os.path.join(OUT, f'z3_{tag}.json'), encoding='utf-8') as f:
        return json.load(f)

say('== pieces (TOOLS sec. S): Z3 h, comm dims (p0,p1), certificate, hull, split')
for tag, x in [('cyc9_013_d7', (1, 0)), ('cyc12_015_d8', (1, 0)), ('cyc8_016_d9', (3, 0)), ('cyc12_014_d9', (1, 0)),
               ('cyc3_012_d8', (2, 1)), ('cyc3_012_d8', (2, 0)), ('cyc5_012_d10', (1, 0)), ('cyc4_013_d10', (1, 2)),
               ('cyc7_012_d8', (2, 0)), ('cyc10_012_d10', (2, 0)), ('cyc10_012_d10', (2, 6)),
               ('cyc7_016_d8', (2, 0)), ('cyc7_016_d8', (6, 0))]:
    v = load(tag)['pieces'][str(x)]
    say(f"  {tag:14s} {x}: h={tuple(v['h'])} delta={v['delta']} comm={v['comm']} certE={v['certE']} certU={v['certU']} [{v['how']}] hull={v['hull']} split={v['split']}")

say('== block details')
for tag, x, y in [('cyc12_014_d9', (1, 0), (1, 4)), ('cyc5_014_d6', (2, 0), (2, 0)), ('cyc7_016_d8', (2, 0), (2, 0)),
                  ('cyc7_016_d8', (2, 0), (2, 1)), ('cyc7_016_d8', (4, 0), (4, 0)), ('cyc9_018_d10', (5, 0), (5, 0))]:
    r = load(tag)
    n, d = r['n'], r['d']
    U = units(lcm(d, n))
    for b in r['blocks']:
        mem = {(gal(t, tuple(b['rep'][0]), d, n), gal(t, tuple(b['rep'][1]), d, n)) for t in U}
        if (x, y) in mem:
            say(f"  {tag} {x}->{y}: rep {tuple(b['rep'][0])}->{tuple(b['rep'][1])} x{b['size']} lb={b['lb']} ub={b['ub']} tool={b['tool']} details={' '.join(b['details'])[:300]}")

say('== Z3 J_K profiles, K = 0..4 (TOOLS sec. S: cyc9_012 d10 (5,3)->(5,6) [57,22,2,0,0]; cyc4_013 d10 (5,1)->(5,3) [241,121,7,0,0])')
for name, d, x, y in [('cyc9_012', 10, (5, 3), (5, 6)), ('cyc4_013', 10, (5, 1), (5, 3))]:
    n, w = parse_family(name)
    for p in pick_primes(d, n):
        pd = PrimeData(name, d, p)
        Px, Py = pd.pieces[x], pd.pieces[y]
        prof = [jtest(Px, Py, [Wa[:k + 1] for Wa in pd.W(x)], [Wa[:k + 1] for Wa in pd.W(y)], k, p) for k in range(5)]
        say(f'  {name} d{d} p={p} {x}->{y}: {prof}  (6 arcs; Z1 used 8)')

say('== closing route: book.txt tool -> Z3 tool (block orbits listed in book.txt)')
cnt = collections.Counter()
for ln in open(os.path.join(OUT, 'z3_compare.txt'), encoding='utf-8'):
    m = re.match(r'\s+OK\s+(\S+)\s+.*book=\[\d+,\d+\] \((.*?)\)\s+Z3=\[\d+,\d+\] \((.*?)\)', ln)
    if m:
        bt = m.group(2).split(' ')[0].split('(')[0]
        cnt[(bt, m.group(3))] += 1
for (bt, zt), c in sorted(cnt.items(), key=lambda q: -q[1]):
    say(f'  {c:4d}  {bt:14s} -> {zt}')
with open(os.path.join(OUT, 'z3_spot.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines) + '\n')
