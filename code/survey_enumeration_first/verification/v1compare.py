# Compare V1 results with the S1 survey outputs (read only AFTER V1 results were saved).
import json, os, re, sys
from math import gcd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTD = os.path.join(BASE, 'outputs')
S1D = os.path.join(os.path.dirname(BASE), 'S1_survey', 'outputs')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v1run import FAMS_A, FAMS_B, FAMS_C, FAMS_D, FAMS_CTRL

lines = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    lines.append(s)
    print(s)


def lcm(a, b):
    return a * b // gcd(a, b)


s1cache = {}


def s1item(label, d):
    if d not in s1cache:
        s1cache[d] = {it['name']: it for it in json.load(open(os.path.join(S1D, f'survey_d{d}_v4.json')))}
    return s1cache[d].get(label)


def lab_s1_to_v1(lab):
    a, chi = map(int, lab.strip('()').split(','))
    return f'C:{chi}' if a == 0 else f'S:{a}:{chi}'


def expand(myJ, T):
    d, n = myJ['d'], myJ['n']
    N = lcm(d, n) if T == 'Y' else n
    units = [k for k in range(1, N) if gcd(k, N) == 1] or [1]

    def act(k, key):
        p = key.split(':')
        if p[0] == 'S':
            return f'S:{(k * int(p[1])) % d}:{(k * int(p[2])) % n}'
        return f'C:{(k * int(p[1])) % n}'
    out = {}
    for b in myJ[T]['blocks']:
        for k in units:
            out[(act(k, b['s']), act(k, b['t']))] = b
    return out


def mykind(b, fixed=False):
    r, lb, ub = (b['fixed'] if fixed else (b['rule'], b['lb'], b['ub']))
    if r == 'R1:NSxNS' or r == 'R1:CMxCM':
        return 'NSxNS/CMxCM'
    if r == 'R4:aut':
        return 'aut'
    if r.startswith('R6') or r == 'R1:pure-res':
        if lb != ub:
            return 'residual'
        return 'aut-IVHS' if lb > 0 else 'zero-IVHS'
    return 'zero'


def s1kind(b):
    k = b.get('kind', '?')
    if k in ('NSxNS', 'CMxCM'):
        return 'NSxNS/CMxCM'
    return k


whys = {}
special = []
fams = []
for f in FAMS_A + FAMS_B + FAMS_C + FAMS_D + FAMS_CTRL:
    if f not in fams:
        fams.append(f)
tot_block_dis = 0
for label, d in fams:
    my = json.load(open(os.path.join(OUTD, f'v1_{label}_d{d}.json')))
    s1 = s1item(label, d)
    tag = 'A' if (label, d) in FAMS_A else ('B' if (label, d) in FAMS_B else ('C' if (label, d) in FAMS_C else ('D' if (label, d) in FAMS_D else 'ctrl')))
    if (label, d) in FAMS_A and (label, d) in FAMS_B:
        tag = 'A+B'
    if s1 is None:
        log(f'== {label} d{d} [{tag}]: not in S1 survey_d{d}_v4.json (control)')
        log(f'   V1 Y={my["Y"]["lb"]}..{my["Y"]["ub"]} C={my["C"]["lb"]}..{my["C"]["ub"]}')
        continue
    log(f'== {label} d{d} [{tag}]')
    # family data
    s1perms = [''.join(map(str, p)) for p in s1.get('normalizer_perms', []) if list(p) != [0, 1, 2]]
    log(f'   |S| V1={my["nS"]} S1={s1["nmon"]}; moduli V1={my["moduli"]} S1={s1.get("moduli_both")}; Gmax V1={my["gmax"]} S1={s1["C"].get("Gmax_order")}; perms V1={my["perms"]} S1={s1perms}')
    # pieces
    bad = []
    s1p = {lab_s1_to_v1(k): tuple(v) for k, v in s1['pieces'].items() if sum(v)}
    myp = {}
    for T in ('Y', 'C'):
        for o in my[T]['orbits']:
            for k, h in o['h'].items():
                myp[k] = tuple(h)
    for k in set(s1p) | set(myp):
        if s1p.get(k) != myp.get(k):
            bad.append((k, s1p.get(k), myp.get(k)))
    log(f'   pieces: {len(myp)} nonzero (V1) vs {len(s1p)} (S1); Hodge-number mismatches: {len(bad)} {bad[:6]}')
    # totals
    log(f'   Y: V1 {my["Y"]["lb"]}..{my["Y"]["ub"]} closed={my["Y"]["closed"]} noIVHS={my["Y"]["closed_without_ivhs"]} | fixed {my["Y"]["fixed"]["lb"]}..{my["Y"]["fixed"]["ub"]} closed={my["Y"]["fixed"]["closed"]} noIVHS={my["Y"]["fixed"]["closed_without_ivhs"]} | S1 {s1["Y_lb"]}..{s1["Y_ub"]} fully_proved_Y={s1["fully_proved_Y"]} by_rules={s1["fully_proved_Y_R"]}')
    log(f'   C: V1 {my["C"]["lb"]}..{my["C"]["ub"]} | fixed {my["C"]["fixed"]["lb"]}..{my["C"]["fixed"]["ub"]} | S1 {s1["C_lb"]}..{s1["C_ub"]} decided={s1.get("C_decided")}')
    log(f'   orbit counts Y: V1 {my["Y"]["orbit_counts"]} S1 {s1["orbits_Y"]}; C: V1 {my["C"]["orbit_counts"]} S1 {s1["orbits_C"]}')
    # orbits
    myorb = {}
    for T in ('Y', 'C'):
        for o in my[T]['orbits']:
            myorb[frozenset(o['key'])] = o
    for o in s1['orbits']:
        mem = frozenset(lab_s1_to_v1(e[0]) for e in o['emb'])
        mo = myorb.get(mem)
        if mo is None:
            log(f'   ORBIT MISMATCH (membership) S1 {sorted(mem)}')
            continue
        k1 = {'residual': 'res'}.get(o['cls'], o['cls']).lower()
        k2 = mo['kind'].lower()
        if k1 != k2:
            att = [(a['member'], a['p'], a['hull'], a.get('l4', {}).get('cert') if a.get('l4') else None) for a in mo.get('attempts', []) if a['member'] == mo.get('cert')][:1]
            log(f'   orbit class differs: {sorted(mem)[:3]}.. S1={o["cls"]} ({o.get("why", "")}) V1={mo["kind"]} cert={mo.get("cert")} {att}')
        elif k1 == 'rigid' and bool(o.get('galois_rigid')) != bool(mo.get('galois_rigid')):
            log(f'   galois_rigid differs: {sorted(mem)[:3]}.. S1={o.get("galois_rigid")} V1={mo.get("galois_rigid")} V1why={mo.get("grigid_why")}')
    # blocks
    for T in ('Y', 'C'):
        mine = expand(my, T)
        cross, dis = {}, []
        n_s1 = 0
        for a, bl in s1['blocks'].items():
            if (T == 'C') != (a == '0'):
                continue
            for b in bl:
                n_s1 += 1
                w = re.sub(r'[-0-9]+', '#', b.get('why', ''))[:90]
                whys[(b.get('kind'), w)] = whys.get((b.get('kind'), w), 0) + 1
                key = (lab_s1_to_v1(b['s']), lab_s1_to_v1(b['t']))
                m = mine.get(key)
                if m is None:
                    dis.append(('missing in V1', key, b.get('kind'), b.get('lb'), b.get('ub')))
                    continue
                if 'Higgs-kernel subspaces differ' in b.get('why', '') or 'both orbits simple' in b.get('why', ''):
                    special.append((label, d, key, b.get('kind'), b.get('why', '')[:60], 'V1 old:', m['rule'], m['lb'], m['ub'], 'V1 fixed:', tuple(m['fixed']), 'V1 zero_rules', m.get('zero_rules'), m.get('zero_rules_fixed')))
                ck = (s1kind(b), mykind(m))
                cross[ck] = cross.get(ck, 0) + 1
                if (b.get('lb'), b.get('ub')) != (m['lb'], m['ub']):
                    dis.append((key, 'S1', b.get('kind'), b.get('lb'), b.get('ub'), b.get('why', '')[:80], 'V1', m['rule'], m['lb'], m['ub']))
        extra = len(mine) - n_s1
        tot_block_dis += len(dis)
        log(f'   {T}-blocks: S1 {n_s1}, V1 {len(mine)}; (lb,ub) disagreements: {len(dis)}; kind cross-tab (S1 kind, V1 kind): {dict(sorted(cross.items()))}')
        for x in dis[:12]:
            log('      DIS', x)
    o5 = my['Y']['only_r5_or_grigid'] + my['C']['only_r5_or_grigid']
    if o5:
        log(f'   V1 blocks closed only by R5/galois_rigid rules: {[(o["s"], o["t"], o["size"], o["zero_rules"], "fixed:", o["zero_rules_fixed"], "ub1:", o["ub1_if_rule_dropped"]) for o in o5]}')
log('')
log('== S1 block kinds and reasons (all compared families), counts:')
for (k, w), c in sorted(whys.items(), key=lambda x: -x[1]):
    log(f'   {c:6d}  {k:12s} {w}')
log(f'total (lb,ub) block disagreements: {tot_block_dis}')
log('== S1 blocks closed by S1 rule 5 (Higgs-kernel subspaces) or rule 3 (both orbits simple), with V1 classification:')
for x in special:
    log('   ', x)
with open(os.path.join(OUTD, 'compare_S1.txt'), 'w') as f:
    f.write('\n'.join(lines) + '\n')
