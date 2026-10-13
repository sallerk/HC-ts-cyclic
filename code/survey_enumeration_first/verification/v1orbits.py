import json, sys, os
OUTD = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'outputs')
for name in sys.argv[1:]:
    J = json.load(open(os.path.join(OUTD, f'v1_{name}.json')))
    print('==', name, 'primes', J['primes'], 'perms', J['perms'], 'moduli', J['moduli'])
    for T in ('Y', 'C'):
        R = J[T]
        print(f' {T}: lb={R["lb"]} ub={R["ub"]} rules={R["rule_counts"]}')
        for o in R['orbits']:
            att = [(a['member'], a['p'], a['hull'], (a['l4']['cert'], a['l4']['minrank_ge'], a['l4']['k']) if a.get('l4') else None) for a in o.get('attempts', [])]
            print(f"   {o['kind']:5s} cert={o.get('cert')} grigid={o.get('galois_rigid')} h={o['h']} att={att}")
        for b in R['blocks']:
            if b['rule'].startswith('R6') or b['rule'] == 'R1:pure-res' or b.get('lb_aut', 1) > 1:
                print('   block', b['s'], '->', b['t'], 'x', b['size'], b['hs'], b['ht'], b['ks'], b['kt'], 'aut', b['aut'], b['rule'], b['lb'], b['ub'], 'zr', b.get('zero_rules'), 'ub1', [v[3] for v in b.get('ub1', [])][:6])
