import sys, re, os, json
sys.stdout.reconfigure(encoding='utf-8')
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'outputs')
cnt = {}
for sub in [os.path.join('survey_second', 'main'), os.path.join('survey_second', 'further_families')]:
    od = os.path.join(BASE, sub)
    for fn in sorted(os.listdir(od)):
        m = re.match(r'z3_(cyc\S+)_d(\d+)\.json$', fn)
        if m:
            r = json.load(open(os.path.join(od, fn), encoding='utf-8'))
            k = len(set(r['primes']))
            cnt[k] = cnt.get(k, 0) + 1
            if k != 2:
                print(sub, fn, r['primes'])
print('distinct-prime counts per run:', cnt)
