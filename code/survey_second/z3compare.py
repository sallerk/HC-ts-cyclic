# z3compare.py -- compare Z3 results with Z1 z1_final_table.txt and Z2b book.txt (run only after Z3 results were saved).
import sys, os, re, json, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf-8')
from z3lib import *
from z3run import FAMILIES, CONTROLS, PrimeData

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'outputs')
T24 = os.path.join(HERE, '..', '..')
BOOK = os.path.join(T24, 'Z2b_tools_writeup', 'outputs', 'book.txt')
Z1T = os.path.join(T24, 'Z1_residual_triage', 'outputs', 'z1_final_table.txt')

def load(name, d):
    with open(os.path.join(OUT, f'z3_{name}_d{d}.json'), encoding='utf-8') as f:
        return json.load(f)

mine = {}
for (name, d) in FAMILIES + CONTROLS:
    r = load(name, d)
    L = lcm(d, r['n']); U = units(L)
    idx = {}
    for b in r['blocks']:
        x, y = tuple(b['rep'][0]), tuple(b['rep'][1])
        for t in U:
            idx[(gal(t, x, d, r['n']), gal(t, y, d, r['n']))] = b
    mine[f'{name}_d{d}'] = (r, idx)

lines = []
def say(s):
    lines.append(s); print(s)

# ---- family Y from book SUMMARY and z1 table
bookY, z1Y = {}, {}
book = open(BOOK, encoding='utf-8').read().splitlines()
for ln in book:
    m = re.match(r'^(cyc\S+_d\d+)\s+Y=(\d+)\.\.(\d+) (CLOSED|OPEN)', ln)
    if m:
        bookY[m.group(1)] = (int(m.group(2)), int(m.group(3)), m.group(4))
z1 = open(Z1T, encoding='utf-8').read().splitlines()
for ln in z1:
    m = re.match(r'^(cyc\S+_d\d+)\s+Y=(\d+)\.\.(\d+)\s+(CLOSED|open)', ln)
    if m:
        z1Y[m.group(1)] = (int(m.group(2)), int(m.group(3)), m.group(4))
say('== per-family Y: Z3 vs Z2b book.txt (final) vs Z1 z1_final_table.txt (before the curve rules)')
nagree = 0
for tag, (r, idx) in mine.items():
    my = (r['Ylb'], r['Yub'])
    bk = bookY.get(tag); zz = z1Y.get(tag)
    ok = bk is not None and (bk[0], bk[1]) == my
    nagree += ok
    zin = zz is not None and zz[0] <= my[0] and my[1] <= zz[1]
    say(f'  {tag:18s} Z3={my[0]}..{my[1]}  book={bk[0]}..{bk[1]} {bk[2]}  agree={ok}   z1table={zz[0]}..{zz[1]} {zz[2]}  Z3 inside z1 range={zin}')
say(f'  agreement with book.txt: {nagree} of {len(mine)}')

# ---- planes
say('== curve planes: book vs Z3')
cur = None
for ln in book:
    m = re.match(r'^(cyc\S+_d\d+)\s+Y=.*planes=(\[.*?\]) notes', ln)
    if m:
        tag = m.group(1)
        bp = eval(m.group(2))
        r = mine[tag][0]
        zp = [(p['i'], tuple(p['M']), p['g']) for p in r['planes']]
        same = [(a, tuple(b), c) for (a, b, c, *_rest) in bp] == zp
        if not same or bp:
            say(f'  {tag:18s} book={bp} Z3={zp} same(i,M,g)={same}')

# ---- block level: every block listed in book.txt
say('== block level: every block listed in book.txt (special-tool closures and OPEN blocks)')
pat = re.compile(r'^\s+(OPEN )?S:(\d+):(\d+)->S:(\d+):(\d+) size=(\d+)\s*(?:h=\((.*?)\)->\((.*?)\) )?(.*)$')
cur = None
nb, ndis = 0, 0
for ln in book:
    m0 = re.match(r'^(cyc\S+_d\d+)\s+Y=', ln)
    if m0:
        cur = m0.group(1); continue
    m = pat.match(ln)
    if not m or cur is None:
        continue
    isopen = bool(m.group(1))
    x = (int(m.group(2)), int(m.group(3))); y = (int(m.group(4)), int(m.group(5)))
    size = int(m.group(6)); rest = m.group(9)
    if isopen:
        mm = re.search(r'\[(\d+),(\d+)\]', rest); blb, bub = int(mm.group(1)), int(mm.group(2)); btool = 'OPEN'
    else:
        mm = re.match(r'-> (\d+)\s*\[(.*?)\]', rest)
        if mm:
            blb = bub = int(mm.group(1)); btool = mm.group(2)
        else:
            mm = re.match(r'closed by (\S+) at (\d+)', rest)
            blb = bub = int(mm.group(2)); btool = mm.group(1)
    r, idx = mine[cur]
    b = idx.get((x, y))
    nb += 1
    if b is None:
        say(f'  {cur} {x}->{y}: NOT FOUND in Z3'); ndis += 1; continue
    hx = r['pieces'][str(x)]['h'] if str(x) in r['pieces'] else None
    hy = r['pieces'][str(y)]['h'] if str(y) in r['pieces'] else None
    hok = True
    if m.group(7):
        hok = (tuple(int(v) for v in m.group(7).split(',')) == tuple(hx) and tuple(int(v) for v in m.group(8).split(',')) == tuple(hy))
    agree = (b['lb'], b['ub']) == (blb, bub) and b['size'] == size and hok
    ndis += (not agree)
    say(f"  {'OK ' if agree else 'DIS'} {cur:18s} {x}->{y} x{size}/{b['size']}  book=[{blb},{bub}] ({btool[:40]})  Z3=[{b['lb']},{b['ub']}] ({(b['tool'] or '').split('@')[0]})  h-match={hok}")
say(f'  blocks compared: {nb}, disagreements: {ndis}')

# ---- z1_final_table open blocks (pre-curve ranges): Z3 value must lie inside
say('== z1_final_table OPEN blocks (pre-curve [lb,ub]) vs Z3 final value')
nz, nout = 0, 0
for ln in z1:
    m = re.match(r'^\s+(cyc\S+_d\d+)\s+S:(\d+):(\d+)->S:(\d+):(\d+) size=(\d+) .*?lb=(\d+) ub=(\d+)', ln)
    if not m:
        continue
    tag = m.group(1)
    x = (int(m.group(2)), int(m.group(3))); y = (int(m.group(4)), int(m.group(5)))
    lb, ub = int(m.group(7)), int(m.group(8))
    r, idx = mine[tag]
    b = idx[(x, y)]
    inside = lb <= b['lb'] and b['ub'] <= ub
    nz += 1; nout += (not inside)
    say(f"  {'in ' if inside else 'OUT'} {tag:18s} {x}->{y} x{m.group(6)}  z1=[{lb},{ub}]  Z3=[{b['lb']},{b['ub']}] {(b['tool'] or '').split('@')[0]} | lbsrc {b['lbsrc']}")
say(f'  z1 open blocks: {nz}, Z3 value outside z1 range: {nout}')

# ---- |H| vs |Gamma image| for all runs
bad = [tag for tag, (r, _) in mine.items() if r['gmax'][0] != r['gmax'][1]]
say(f'== |H| = |Gamma image| in all 33 runs: {not bad} {bad}')

# ---- own J dims by order at cyc8_013 d6 (1,2)->(1,4) and (1,6)->(1,4) (book prints [7,5,3,1,0] at 50329/50377)
say('== Z3 J_K, K = 0..4, at cyc8_013 d6 (book: [7, 5, 3, 1, 0] at its primes)')
for p in pick_primes(6, 8):
    pd = PrimeData('cyc8_013', 6, p)
    for (x, y) in [((1, 2), (1, 4)), ((1, 6), (1, 4)), ((1, 4), (1, 2))]:
        Px, Py = pd.pieces[x], pd.pieces[y]
        Wx, Wy = pd.W(x), pd.W(y)
        prof = [jtest(Px, Py, [Wa[:k + 1] for Wa in Wx], [Wa[:k + 1] for Wa in Wy], k, p) for k in range(5)]
        say(f'  p={p} {x}->{y}: {prof}')

with open(os.path.join(OUT, 'z3_compare.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines) + '\n')
