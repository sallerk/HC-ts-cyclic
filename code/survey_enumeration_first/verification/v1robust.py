# Compare the main V1 run (primes from 50000) with a second run at fresh primes/points (V1_START=70000),
# and print the galois_rigid-only blocks of cyc12_016 d7.
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v1run import FAMS_A, FAMS_B, FAMS_C, FAMS_D, FAMS_CTRL
OUTD = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'outputs')
lines = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    lines.append(s)
    print(s)


fams = []
for f in FAMS_A + FAMS_B + FAMS_C + FAMS_D + FAMS_CTRL:
    if f not in fams:
        fams.append(f)
ndiff = 0
for label, d in fams:
    A = json.load(open(os.path.join(OUTD, f'v1_{label}_d{d}.json')))
    B = json.load(open(os.path.join(OUTD, f'v1_{label}_d{d}_s70000.json')))
    diffs = []
    for T in ('Y', 'C'):
        for key in ('lb', 'ub', 'closed', 'closed_without_ivhs', 'orbit_counts'):
            if A[T][key] != B[T][key]:
                diffs.append((T, key, A[T][key], B[T][key]))
        for key in ('lb', 'ub', 'closed', 'closed_without_ivhs'):
            if A[T]['fixed'][key] != B[T]['fixed'][key]:
                diffs.append((T, 'fixed.' + key, A[T]['fixed'][key], B[T]['fixed'][key]))
        ka = {(b['s'], b['t']): (b['rule'], b['lb'], b['ub'], tuple(b['fixed'])) for b in A[T]['blocks']}
        kb = {(b['s'], b['t']): (b['rule'], b['lb'], b['ub'], tuple(b['fixed'])) for b in B[T]['blocks']}
        nb = sum(1 for k in ka if ka[k] != kb.get(k))
        if nb:
            diffs.append((T, 'blocks differing (rule/lb/ub)', nb))
    ndiff += len(diffs)
    log(f'{label:14s} d{d} primes {A["primes"]} vs {B["primes"]}: ' + ('identical' if not diffs else str(diffs)))
log(f'families with differences: {ndiff}')
J = json.load(open(os.path.join(OUTD, 'v1_cyc12_016_d7.json')))
log('cyc12_016 d7: blocks closed only by the galois_rigid rule (R3/R3b), with the first-order bound if the rule is dropped:')
for o in J['Y']['only_r5_or_grigid']:
    log('   ', o)
for o in J['Y']['orbits']:
    if o['kind'] == 'rigid':
        log('   orbit cert', o['cert'], 'galois_rigid', o.get('galois_rigid'), 'reasons', o.get('grigid_why'))
with open(os.path.join(OUTD, 'robustness_s70000.txt'), 'w') as f:
    f.write('\n'.join(lines) + '\n')
