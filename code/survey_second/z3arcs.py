# z3arcs.py -- J_K profile vs number of arcs (diagnosis of the cyc4_013 d10 (5,1)->(5,3) profile difference).
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf-8')
from z3lib import *
import z3run

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'outputs')
lines = []
for na in (6, 8, 10):
    z3run.NARCS = na
    for name, d, x, y in [('cyc4_013', 10, (5, 1), (5, 3)), ('cyc9_012', 10, (5, 3), (5, 6))]:
        n, w = parse_family(name)
        for p in pick_primes(d, n):
            pd = z3run.PrimeData(name, d, p)
            prof = [jtest(pd.pieces[x], pd.pieces[y], [Wa[:k + 1] for Wa in pd.W(x)], [Wa[:k + 1] for Wa in pd.W(y)], k, p)
                    for k in range(5)]
            s = f'arcs={na:2d} {name} d{d} p={p} {x}->{y}: J_0..J_4 = {prof}  (h_x={pd.pieces[x].h}, h_y={pd.pieces[y].h})'
            print(s, flush=True); lines.append(s)
with open(os.path.join(OUT, 'z3_arcs.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines) + '\n')
