"""Summary tables for the S1 survey JSONs.  Usage: python -I -B summarize.py TAG d1 d2 ...  (reads survey_d<d><TAG>.json)"""
import json
import os
import sys
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, '..', 'outputs'))
tag = sys.argv[1] if sys.argv[1] != '-' else ''
ds = [int(x) for x in sys.argv[2:]]
L = []
P = print


def decomposable(r):
    """F = c x_i^d + g(other two variables) for some i (x_i occurs only as x_i^d)."""
    d, n, w = r['d'], r['n'], r['w']
    S = [(a, b, d - a - b) for a in range(d, -1, -1) for b in range(d - a, -1, -1)
         if (w[1] * b + w[2] * (d - a - b)) % n == 0]
    return [i for i in range(3) if all(e[i] in (0, d) for e in S)]


def recompute(r):
    """ub with R3 blocks capped at 1 (rigid source, bijective block orbit: a nonzero Hodge map is an isomorphism and
    phi^-1 phi' lies in the K-linear Hodge endomorphisms = K, so the block has dim <= 1)."""
    out = {}
    for part, sel in (('Y', lambda a: a != '0'), ('C', lambda a: a == '0')):
        tot, cnt = 0, Counter()
        exc = Counter()
        for a, bl in r['blocks'].items():
            if not sel(a):
                continue
            for b in bl:
                u = b['ub']
                if b['kind'] == 'residual':
                    key = ('R3' if b.get('why', '').startswith('R3') else
                           'aut-nonrigid' if b.get('why', '').startswith('automorphism') else
                           'incomparable' if 'incomparable' in b.get('why', '') else 'resorb')
                    if key == 'R3':
                        u = 1 if u is None else min(u, 1)
                    cnt[key] += 1
                    exc[key] += (u - b['lb']) if u is not None else 0
                if tot is not None:
                    tot = None if u is None else tot + u
        out[part] = (tot, cnt, exc)
    return out


def rng(lb, ub):
    return str(lb) if lb == ub else "%s..%s" % (lb, ub)


proved_R, proved_I, notp, shaky = [], [], [], []
for d in ds:
    fn = os.path.join(OUT, "survey_d%d%s.json" % (d, tag))
    if not os.path.exists(fn):
        P("missing", fn)
        continue
    allr = json.load(open(fn, encoding='utf-8'))
    ok = [r for r in allr if r.get('status') == 'ok']
    drop = Counter(r['status'] for r in allr if r.get('status') != 'ok')
    P("== d=%d: %d monomial sets (n<=%d), kept %d (smooth, moduli>0); dropped: %s" % (d, len(allr), 2 * d, len(ok),
                                                                                   dict(drop)))
    P("family          n  w        |S| mod |Gmax| perms   Y:rig/pure/res C:rig/CM/res  R1-R3  +IVHS  Y-dim        "
      "C-dim     resY")
    for r in ok:
        oy, oc = r['orbits_Y'], r['orbits_C']
        perms = [''.join(map(str, p)) for p in r['normalizer_perms'] if p != [0, 1, 2]]
        nresY = sum(1 for b in r['residual_blocks'] if b['a'] != 0)
        r['dec'] = decomposable(r)
        rc = recompute(r)
        r['Y_ub_raw'], r['C_ub_raw'] = r['Y_ub'], r['C_ub']
        r['Y_ub'], r['C_ub'] = rc['Y'][0], rc['C'][0]
        r['res_counts'] = rc
        P("%-15s%s %2d %-8s %3d %3d %4d %-7s %2d/%2d/%2d       %2d/%2d/%2d      %-5s  %-5s  %-12s %-9s %d" % (
            r['name'], 'D' if r['dec'] else ' ', r['n'], ''.join(map(str, r['w'])) if max(r['w']) < 10 else '-'.join(map(str, r['w'])),
            r['nmon'], r['moduli'], r['C']['Gmax_order'], ','.join(perms) or '-',
            oy.get('rigid', 0), oy.get('pure', 0), oy.get('residual', 0),
            oc.get('rigid', 0), oc.get('CM', 0), oc.get('residual', 0),
            'YES' if r['fully_proved_Y_R'] else 'no', 'YES' if r['fully_proved_Y'] else 'no',
            rng(r['Y_lb'], r['Y_ub']), rng(r['C_lb'], r['C_ub']) + ('' if r['C_decided'] else '*'), nresY))
        if r['fully_proved_Y_R']:
            proved_R.append(r)
        elif r['fully_proved_Y']:
            proved_I.append(r)
        else:
            notp.append(r)
        fl = []
        for k in ('prime_disagree_rigid', 'prime_disagree_ub1', 'inconsistent', 'stab_span_notes'):
            if r.get(k):
                fl.append("%s=%s" % (k, r[k][:3]))
        if r['hilbert_mismatch']:
            fl.append('hilbert_mismatch=%d' % r['hilbert_mismatch'])
        for a, bl in r['blocks'].items():
            for b in bl:
                if b.get('flag'):
                    fl.append('a=%s %s->%s %s' % (a, b['s'], b['t'], b['flag']))
                if b['kind'] == 'zero' and b.get('why', '').startswith('R3'):
                    if not all(x['rank_maximal'] for x in b['vhs']):
                        fl.append('a=%s %s->%s R3-zero with non-maximal mod-p ranks (generic-rank assumption)' % (
                            a, b['s'], b['t']))
                if 'incomparable' in b.get('why', ''):
                    fl.append('a=%s %s->%s stabilizers incomparable' % (a, b['s'], b['t']))
        if fl:
            shaky.append((r['name'], d, fl))


def hdg(r):
    C = rng(r['C_lb'], r['C_ub']) + ('' if r['C_decided'] else ' [C-part ub first order, inf]')
    tot = (1 + r['Y_lb'] + r['C_lb'], 1 + r['Y_ub'] + (r['C_ub'] if r['C_ub'] is not None else 0))
    return "Y=%s C=%s  dim Hdg^{2,2}(X_FF) = 1+Y+C = %s" % (rng(r['Y_lb'], r['Y_ub']), C, rng(*tot))


P("\n== FULLY PROVED by R1-R3 alone (no residual block in any j of the Y-part)")
for r in proved_R:
    P("  %-15s d=%d moduli=%d  %s" % (r['name'], r['d'], r['moduli'], hdg(r)))
P("\n== FULLY PROVED only after first-order (IVHS) closure of the R4 blocks (ub1 = lb at 2 primes)")
for r in proved_I:
    nI = sum(1 for bl in r['blocks'].values() for b in bl if b['kind'] in ('aut-IVHS', 'zero-IVHS')
             and b['s'] and int(r['blocks'] and 1))
    kinds = Counter(b['kind'] for a, bl in r['blocks'].items() if a != '0' for b in bl)
    P("  %-15s d=%d moduli=%d  %s  (Y closures: aut-IVHS %d, zero-IVHS %d)" % (
        r['name'], r['d'], r['moduli'], hdg(r), kinds.get('aut-IVHS', 0), kinds.get('zero-IVHS', 0)))
P("\n== NOT proved: residual Y-blocks (first-order bound ub1 = Galois-min over the block orbit, 2 primes)")
for r in notp:
    P("  %-15s%s d=%d moduli=%d Y-dim %s | residual Y-blocks by reason %s, sum(ub-lb) %s" % (
        r['name'], ' D' if r['dec'] else '', r['d'], r['moduli'], rng(r['Y_lb'], r['Y_ub']),
        dict(r['res_counts']['Y'][1]), dict(r['res_counts']['Y'][2])))
    for b in r['residual_blocks']:
        if b['a'] == 0:
            continue
        P("     a=%d %s%s -> %s%s %s lb=%s ub1=%s : %s" % (b['a'], b['s'], b['hs'], b['t'], b['ht'], b['cls'], b['lb'],
                                                       b['ub1'], b['why']))
P("\n== flags (possible shaky points)")
for nm, d, fl in shaky:
    P("  %s d=%d: %s" % (nm, d, " | ".join(fl[:6]) + (" ..." if len(fl) > 6 else "")))
