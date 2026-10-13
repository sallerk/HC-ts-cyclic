"""Compare S1 survey totals with B's famk chain numbers (read-only B JSON), and optionally list blocks.
Usage: python -I -B compareB.py SURVEY.json FAMNAME [BJSON] [show=a1,a2]
chain r of famk <-> x-piece i = (-4-r) mod d <-> u-character a = i+1 mod d (a = 0: invariant chain = H^1(C)xH^1(C))."""
import json
import sys

sys.stdout.reconfigure(encoding='utf-8')
surv = json.load(open(sys.argv[1], encoding='utf-8'))
name = sys.argv[2]
bj = None
show = []
for a in sys.argv[3:]:
    if a.startswith('show='):
        show = [int(x) for x in a[5:].split(',')]
    else:
        bj = json.load(open(a, encoding='utf-8'))
res = [r for r in surv if r['name'] == name][0]
d = res['d']
print("%s d=%d status=%s fully_proved_R=%s fully_proved_IVHS=%s" % (name, d, res['status'], res.get('fully_proved_Y_R'),
                                                                  res.get('fully_proved_Y')))
if bj:
    print("B file: family=%s p=%s seed=%s" % (bj['family'], bj['p'], bj['seed']))
    print(" r  i  a | B: unknowns K Kgal LB | S1: lb ub (ub1 galmin-sum raw-sum if computed)")
    for r in range(d):
        c = bj['chains'][str(r)]
        i = (-4 - r) % d
        a = (i + 1) % d
        t = res['totals'].get(str(a), {})
        print("%2d %2d %2d | %6d %5d %5d %5d | %s %s %s %s" % (r, i, a, c['unknowns'], c['K'], c['Kgal'], c['LB'],
                                                          t.get('lb'), t.get('ub'), t.get('ub1_galmin_sum', ''),
                                                          t.get('ub1_raw_sum', '')))
for a in show:
    print("blocks a=%d:" % a)
    for b in res['blocks'][str(a)]:
        print("  %s%s -> %s%s cls=%s kind=%s lb=%s ub=%s %s %s %s" % (
            b['s'], b['hs'], b['t'], b['ht'], b['cls'], b['kind'], b['lb'], b['ub'], b.get('why', ''),
            b.get('perms', ''), b.get('flag', '')))
for o in res['orbits']:
    print(" orbit %d w%d %s %s %s" % (o['id'], o['weight'], o['cls'], o['emb'],
                                      [(t['emb'], t['test'], t.get('ranks', t.get('res')), t.get('hull_by_order'),
                                        t['pass']) for t in o.get('tried', [])]))
print("stab span notes:", res.get('stab_span_notes'))
