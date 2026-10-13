import json, sys, os
OUTD = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'outputs')
for name in sys.argv[1:]:
    J = json.load(open(os.path.join(OUTD, f'v1_{name}.json')))
    print('==', name, 'primes', J['primes'], 'perms', J['perms'])
    for T in ('Y', 'C'):
        R = J[T]
        print(f' {T}: lb={R["lb"]} ub={R["ub"]} rules={R["rule_counts"]}')
        for o in R['orbits']:
            if o['kind'] in ('res',) or '-v' in sys.argv[0]:
                print('   res orbit', o['key'], o['h'], [(a['member'], a['p'], a['hull'], a.get('l4')) for a in o.get('attempts', [])])
            if o['kind'] == 'rigid' and not o.get('galois_rigid'):
                print('   rigid not galois-rigid', o['key'], o.get('grigid_why'))
        for b in R['blocks']:
            if not b['closed'] or b['rule'].startswith('R6') or b['rule'] == 'R1:pure-res':
                print('   block', b['s'], '->', b['t'], 'x', b['size'], b['hs'], b['ht'], b['ks'], b['kt'], 'aut', b['aut'], b['rule'], b['lb'], b['ub'], 'ub1', b.get('ub1'), 'r5n', b.get('r5naive'))
        for c in R['conflicts']:
            print('   CONFLICT', c)
