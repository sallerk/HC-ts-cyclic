# z3compare2.py -- every block orbit: Z3 vs Z1 z1o3_<family>.json (V1/ord3 stage) overridden by Z2b book.txt (final).
import sys, os, re, json, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf-8')
from z3lib import lcm, units, gal
from z3run import FAMILIES, CONTROLS

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'outputs')
T24 = os.path.join(HERE, '..', '..')
Z1D = os.path.join(T24, 'Z1_residual_triage', 'outputs')
BOOK = os.path.join(T24, 'Z2b_tools_writeup', 'outputs', 'book.txt')

def lab(s):
    _, a, c = s.split(':')
    return (int(a), int(c))

def find_blocks(o, acc, path=''):
    if isinstance(o, dict):
        if {'s', 't', 'size', 'lb', 'ub'} <= set(o):
            acc.append((path, o))
        for k, v in o.items():
            find_blocks(v, acc, path + '/' + str(k))
    elif isinstance(o, list):
        for v in o:
            find_blocks(v, acc, path)
    return acc

# book values
book = {}
cur = None
pat = re.compile(r'^\s+(OPEN )?S:(\d+):(\d+)->S:(\d+):(\d+) size=(\d+)\s*(?:h=\((.*?)\)->\((.*?)\) )?(.*)$')
for ln in open(BOOK, encoding='utf-8').read().splitlines():
    m0 = re.match(r'^(cyc\S+_d\d+)\s+Y=', ln)
    if m0:
        cur = m0.group(1); continue
    m = pat.match(ln)
    if not m or cur is None:
        continue
    x = (int(m.group(2)), int(m.group(3))); y = (int(m.group(4)), int(m.group(5)))
    rest = m.group(9)
    if m.group(1):
        mm = re.search(r'\[(\d+),(\d+)\]', rest); v = (int(mm.group(1)), int(mm.group(2)))
    else:
        mm = re.match(r'-> (\d+)', rest)
        v = (int(mm.group(1)),) * 2 if mm else (int(re.match(r'closed by \S+ at (\d+)', rest).group(1)),) * 2
    book[(cur, x, y)] = v

lines = []
def say(s):
    lines.append(s); print(s)

tot_cmp, tot_dis, tot_h, tot_missing = 0, 0, 0, 0
rulestat = collections.Counter()
for (name, d) in FAMILIES + CONTROLS:
    tag = f'{name}_d{d}'
    with open(os.path.join(OUT, f'z3_{tag}.json'), encoding='utf-8') as f:
        r = json.load(f)
    n = r['n']; U = units(lcm(d, n))
    myidx = {}
    for b in r['blocks']:
        x, y = tuple(b['rep'][0]), tuple(b['rep'][1])
        for t in U:
            myidx[(gal(t, x, d, n), gal(t, y, d, n))] = b
    fn = os.path.join(Z1D, f'z1o3_{tag}.json')
    if not os.path.exists(fn):
        fn = os.path.join(Z1D, f'z1_{tag}.json')
    with open(fn, encoding='utf-8') as f:
        z = json.load(f)
    entries = find_blocks(z, [])
    zval = {}
    for path, o in entries:
        key = (lab(o['s']), lab(o['t']))
        if 'open_blocks' in path and key in zval:
            continue
        zval[key] = (o, path)
    # orbit-level: map every Z1 entry to my block orbit
    seen_my = set()
    dis, cmpd, hbad = [], 0, 0
    for key, (o, path) in zval.items():
        b = myidx.get(key)
        if b is None:
            dis.append(f'{key} not a Z3 block'); continue
        bid = (tuple(b['rep'][0]), tuple(b['rep'][1]))
        if bid in seen_my:
            continue
        seen_my.add(bid)
        fin = book.get((tag, key[0], key[1]))
        src = 'book'
        if fin is None:
            fin, src = (o['lb'], o['ub']), 'z1o3:' + str(o.get('rule'))
        rulestat[src.split(':')[0] if src == 'book' else 'z1o3']
        cmpd += 1
        hs, ht = o.get('hs'), o.get('ht')
        if hs is not None and (tuple(hs) != tuple(r['pieces'][str(key[0])]['h']) or tuple(ht) != tuple(r['pieces'][str(key[1])]['h'])):
            hbad += 1
        if (b['lb'], b['ub']) != tuple(fin) or b['size'] != o['size']:
            dis.append(f"{key[0]}->{key[1]} x{o['size']}/{b['size']} Z1/Z2b={list(fin)} [{src}] Z3=[{b['lb']},{b['ub']}] [{(b['tool'] or '').split('@')[0]}]")
    missing = [ (tuple(b['rep'][0]), tuple(b['rep'][1])) for b in r['blocks'] if (tuple(b['rep'][0]), tuple(b['rep'][1])) not in seen_my]
    tot_cmp += cmpd; tot_dis += len(dis); tot_h += hbad; tot_missing += len(missing)
    say(f'{tag:18s} Z3 block orbits {len(r["blocks"]):4d}  matched {cmpd:4d}  Z3-only {len(missing):3d}  disagreements {len(dis)}  h-mismatch {hbad}  ({os.path.basename(fn)})')
    for s in dis[:20]:
        say('    ' + s)
    if missing:
        say('    Z3 block orbits without a Z1 entry: ' + ' '.join(f'{a}->{b}' for a, b in missing[:12]) + (' ...' if len(missing) > 12 else ''))
say(f'TOTAL matched block orbits {tot_cmp}, disagreements {tot_dis}, Hodge-number mismatches {tot_h}, Z3 orbits without Z1 entry {tot_missing}')
with open(os.path.join(OUT, 'z3_compare_allblocks.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines) + '\n')
