import json, sys
fn, fam = sys.argv[1], sys.argv[2]
J = json.load(open(fn))
for item in J:
    if item.get('name') == fam:
        print('keys', list(item.keys()))
        for k in ('normalizer_perms', 'stab_span_notes', 'orbits_Y', 'orbits_C', 'fully_proved_Y', 'fully_proved_Y_R', 'Y_lb', 'Y_ub', 'C_lb', 'C_ub', 'C_decided', 'moduli_both', 'primes', 'C'):
            print(k, '=', json.dumps(item.get(k))[:400])
        for a, v in item['totals'].items():
            print('totals', a, json.dumps(v)[:300])
        for o in item['orbits']:
            print('orbit', json.dumps(o)[:700])
        extra = [k for k in item.keys() if k not in ('pieces', 'orbits', 'totals')]
        print('pieces sample', list(item['pieces'].items())[:6])
        for k in item.keys():
            if k not in ('pieces', 'orbits', 'totals', 'normalizer_perms', 'stab_span_notes', 'orbits_Y', 'orbits_C', 'C'):
                v = item[k]
                if isinstance(v, (list, dict)) and len(json.dumps(v)) > 100:
                    print('BIG', k, type(v).__name__, len(v), json.dumps(v)[:600])
        break
