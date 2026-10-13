import sys, re, os, json
sys.stdout.reconfigure(encoding='utf-8')
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'outputs')
S1 = os.path.join(BASE, 'survey_enumeration_first', 'enumeration', 'summary_v4.txt')
OUT = sys.argv[1] if len(sys.argv) > 1 else None

z = {}
primes = set()
for sub in [os.path.join('survey_second', 'main'), os.path.join('survey_second', 'further_families')]:
    od = os.path.join(BASE, sub)
    for fn in os.listdir(od):
        m = re.match(r'z3_(cyc\S+)_d(\d+)\.json$', fn)
        if m:
            r = json.load(open(os.path.join(od, fn), encoding='utf-8'))
            z[(m.group(1), int(m.group(2)))] = r
            primes.update(r['primes'])

rows = []
d = None
for line in open(S1, encoding='utf-8'):
    m = re.match(r'== d=(\d+):', line)
    if m:
        d = int(m.group(1)); continue
    if line.startswith('== '):
        d = None; continue
    if d is None:
        continue
    m = re.match(r'(cyc\S+)\s+(D\s+)?(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s', line)
    if not m:
        continue
    name, dec, n, w, nmon, mod = m.group(1), bool(m.group(2)), int(m.group(3)), m.group(4), int(m.group(5)), int(m.group(6))
    ws = w.split('-') if '-' in w else list(w)
    rows.append(dict(name=name, dec=dec, d=d, n=n, ws=ws, nmon=nmon, mod=mod))

W = {('cyc12_017', 6): 1017, ('cyc14_018', 7): 288}
bad = []
missing = []
cnt = {}
for r in rows:
    key = (r['name'], r['d'])
    if key in W:
        r['route'] = 'W'; r['Y'] = W[key]
        zz = z.get(key)
        if zz is None or zz['Ylb'] > W[key] or zz['Yub'] != W[key]:
            bad.append(key)
    else:
        zz = z.get(key)
        if zz is None:
            missing.append(key); r['Y'] = None
        else:
            if not zz['closed'] or zz['Ylb'] != zz['Yub']:
                bad.append(key)
            r['Y'] = zz['Ylb']
        if r['n'] == 1:
            r['route'] = 'P'
            if r['Y'] is not None and r['Y'] != r['d'] - 1:
                bad.append(key)
        elif r['ws'] == ['0', '0', '1']:
            r['route'] = 'A'
        else:
            r['route'] = 'C'
    cnt[r['route']] = cnt.get(r['route'], 0) + 1
print('rows', len(rows), cnt)
print('missing', missing)
print('bad', bad)
print('primes', min(primes), max(primes), len(primes))
lines = []
for r in rows:
    Y = '?' if r['Y'] is None else str(r['Y'])
    lines.append('%d & %d & %s & %d & %d & %s & %s%s' % (r['d'], r['n'], ','.join(r['ws']), r['nmon'], r['mod'], Y,
                                                        r['route'], '$^*$' if r['dec'] else ''))
half = (len(lines) + 1) // 2
L, R = lines[:half], lines[half:]
out = []
for i in range(half):
    right = R[i] if i < len(R) else '&&&&&&'
    out.append(L[i] + ' && ' + right + ' \\\\')
txt = '\n'.join(out)
print(txt)
if OUT:
    open(OUT, 'w', encoding='utf-8').write(txt + '\n')
