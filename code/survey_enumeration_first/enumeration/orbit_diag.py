"""Why residual orbits are residual (S1 survey JSONs).  Usage: python -I -B orbit_diag.py TAG d1 d2 ..."""
import json
import os
import sys
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')
OUT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'outputs'))
tag = sys.argv[1]
tot = Counter()
for d in [int(x) for x in sys.argv[2:]]:
    allr = json.load(open(os.path.join(OUT, 'survey_d%d%s.json' % (d, tag)), encoding='utf-8'))
    for r in allr:
        if r.get('status') != 'ok':
            continue
        for o in r['orbits']:
            if o['cls'] != 'residual':
                continue
            tr = o.get('tried', [])
            if not tr:
                why = 'no embedding with h_top in {1,2}' + (' and h02=0' if o['weight'] == 2 else '')
            else:
                kinds = []
                for t in tr:
                    if t['test'] == 'hull1' and 'hull_by_order' not in t:
                        a = int(t['emb'][1:t['emb'].index(',')])
                        kinds.append('hull1 deficient %s, higher order unavailable (w%d, j=%d)' % (
                            t['ranks'], o['weight'], a - 1))
                    elif t['test'] == 'hull1':
                        kinds.append('hull<=3 deficient %s of %s' % (t['hull_by_order'][0], t.get('dim')))
                    else:
                        kinds.append('nondeg2 fails (%s)' % t['res'][0][1])
                why = '; '.join(kinds)
            key = 'w%d: %s' % (o['weight'], why.split(' [')[0].split(' (')[0].split(' %')[0])
            cat = ('no-candidate' if not tr else 'higher-order-unavailable' if 'unavailable' in why and
                   all('unavailable' in k or 'nondeg2' in k for k in why.split('; ')) else
                   'hull3-deficient' if 'hull<=3' in why else 'nondeg2-fails')
            tot[(d, o['weight'], cat)] += 1
            if cat == 'higher-order-unavailable' and o['weight'] == 2:
                print("d=%d %-14s orbit %d %s : %s" % (d, r['name'], o['id'], o['emb'][:3], why))
print("counts (d, weight, category):")
for k in sorted(tot):
    print("  ", k, tot[k])
