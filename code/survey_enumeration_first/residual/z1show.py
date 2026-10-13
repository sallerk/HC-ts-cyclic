# print open blocks / split details / residual orbits of z1_<fam>_d<d>.json files
import os, sys, json
sys.stdout.reconfigure(encoding='utf-8')
OUTD = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'outputs')
mode = sys.argv[1]
for tag in sys.argv[2:]:
    fn = os.path.join(OUTD, f'z1_{tag}.json')
    D = json.load(open(fn))
    print(f"== {tag} Y={D['Ylb']}..{D['Yub']} closed={D['closed']} perms={D['perms']} S={D['S']}")
    if mode in ('orb', 'all'):
        for o in D['orbits']:
            if o['kind'] == 'pure':
                continue
            print(f"  orbit {o['kind']:5s} h={list(o['h'].values())[0]} n={len(o['members'])} certs={o.get('member_certs')} "
                  f"cdims={o.get('cdims')} hull3={o.get('hull3')} hull5={o.get('hull5')} grigid={o.get('galois_rigid')} {o.get('grigid_why')}")
    if mode in ('open', 'all', 'split'):
        for r in D['blocks']:
            if mode == 'split' and 'split' not in r and 'split_error' not in r:
                continue
            if mode != 'split' and not r.get('open'):
                continue
            print(f"  block {r['s']}->{r['t']} size={r['size']} h={r['hs']}->{r['ht']} kinds={r['ks']},{r['kt']} aut={r['aut']} "
                  f"fixed={r['fixed']} z={r['z']} open={r.get('open')} ub1={[v[3] for v in r.get('ub1', [])][:4]}")
            if 'split_error' in r:
                print('   SPLIT ERROR', r['split_error'][-400:])
            if 'split' in r and mode == 'split':
                for dd in (r['split']['detail'] if os.environ.get('ALLD') else r['split']['detail'][:2]):
                    print(f"     {dd['x']}->{dd['y']} p={dd['p']} lb={dd['lb']} ub={dd['ub']} pi={dd['pi_x']},{dd['pi_y']} checks={dd['checks']}")
                    for pp in dd['parts']:
                        print('        ', pp)
