# J-test controls: (1) agreement J vs ord3 on blocks where both were computed; (2) positive controls: J on identity
# blocks (s,s) of the pieces involved in J-closed blocks (must stay >= 1 at every order) and on automorphism blocks.
import os, sys, json, glob, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v1lib import *
from v1run import Fam, kstr, fam_from_label
import jtest
sys.stdout.reconfigure(encoding='utf-8')
OUTD = os.path.join(os.path.dirname(HERE), 'outputs')
parse = lambda s: (s.split(':')[0], int(s.split(':')[1]), int(s.split(':')[2]))
out = []
# (1) agreement
agree = dis = 0
for fn in glob.glob(os.path.join(OUTD, 'jt_*.json')):
    J = json.load(open(fn))
    tag = f"{J['name']}_d{J['d']}"
    Z = json.load(open(os.path.join(OUTD, f'z1o3_{tag}.json')))
    o3 = {(r['s'], r['t']): (r.get('ord3') or {}).get('best') for r in Z['blocks'] if r.get('open')}
    for b in J['blocks']:
        ob = o3.get((b['s'], b['t']))
        if ob is None or b['jbest'] is None:
            continue
        if (ob == 0) == (b['jbest'] == 0):
            agree += 1
        else:
            dis += 1
            out.append(f"  DISAGREE {tag} {b['s']}->{b['t']} ord3={ob} J={b['jbest']}")
out.append(f'(1) J vs ord3 zero/nonzero agreement on block orbits: agree={agree} disagree={dis}')
# (2) positive controls
ctrl = {'cyc3_012_d8': ['S:4:1'], 'cyc4_012_d8': ['S:4:1'], 'cyc4_013_d10': ['S:5:1'], 'cyc5_012_d10': ['S:5:1', 'S:5:2'],
        'cyc6_013_d10': ['S:5:1'], 'cyc9_012_d10': ['S:5:3'], 'cyc7_016_d9': ['S:3:1'], 'cyc8_017_d10': ['S:5:1', 'S:5:2'],
        'cyc8_014_d8': ['S:4:1', 'S:4:2'], 'cyc6_015_d8': ['S:4:1']}
for tag, keys in ctrl.items():
    label, d = tag.rsplit('_d', 1)[0], int(tag.rsplit('_d', 1)[1])
    n, w = fam_from_label(label, d)
    F = Fam(label, d, w, n, start=50000)
    p = F.primes[0]
    for ks in keys:
        x = parse(ks)
        if F.h(x)[0] == 0:
            x = F.act(-1, x)
        rng = random.Random(f'Z1-Jctrl|{tag}|{ks}')
        arcs = [{m: rng.randrange(p) for m in F.S} for _ in range(jtest.NARCS)]
        dims = jtest.jdims(F, p, x, x, arcs)
        out.append(f'(2) identity block {tag} {kstr(x)} h={F.h(x)} J dims by order = {dims}  (expect >= 1)')
# automorphism blocks (both rigid) of cyc5_014 d6 and cyc7_016 d8
for tag in ('cyc5_014_d6', 'cyc7_016_d8'):
    label, d = tag.rsplit('_d', 1)[0], int(tag.rsplit('_d', 1)[1])
    n, w = fam_from_label(label, d)
    F = Fam(label, d, w, n, start=50000)
    p = F.primes[0]
    done = 0
    for x in sorted(F.pieces):
        if x[0] != 'S' or F.h(x)[0] < 1:
            continue
        for y in F.aorbit(x):
            if y != x and F.h(y)[0] == F.h(x)[0] and done < 3:
                rng = random.Random(f'Z1-Jaut|{tag}|{kstr(x)}')
                arcs = [{m: rng.randrange(p) for m in F.S} for _ in range(jtest.NARCS)]
                out.append(f'(2) automorphism block {tag} {kstr(x)}->{kstr(y)} h={F.h(x)} J dims = {jtest.jdims(F, p, x, y, arcs)}  (expect >= 1)')
                done += 1
open(os.path.join(OUTD, 'jcontrol.txt'), 'w', encoding='utf-8').write('\n'.join(out) + '\n')
print('\n'.join(out))
