# z3tests.py -- unit tests for z3lib (agent Z3). Run: python -I -B z3tests.py
import sys, os, random, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
from z3lib import *

OUT = []
def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)

def binom_mod(x, k, p):
    num, den = 1, 1
    for j in range(k):
        num = num * ((x - j) % p) % p
        den = den * (j + 1) % p
    return num * pow(den, p - 2, p) % p

def c_n(D, g, n):
    r = g % n
    if r > D or D < 0:
        return 0
    K = (D - r) // n
    return (K + 1) * (D - r + 1) - n * K * (K + 1) // 2

def make_ring(name, d, p, salt='T'):
    n, w = parse_family(name)
    S = family_S(d, n, w)
    for tries in range(20):
        rng = random.Random(f'Z3|{salt}|{name}|{d}|{p}|{tries}')
        F = {m: rng.randrange(1, p) for m in S}
        R = Ring(d, n, w, F, p)
        if not R.koszul_check():
            return R, rng, tries
    raise RuntimeError('no smooth point')

def pieces_of(R):
    P = {}
    for a in range(1, R.d):
        for chi in range(R.n):
            x = Piece(R, a, chi)
            if x.dim:
                P[(a, chi)] = x
    return P

fails = 0
p = pick_primes(8, 8)[0]
say('test prime', p, 'primes for (d,n)=(10,20):', pick_primes(10, 20))

# 1. Hilbert function = Koszul; total dimension; THEOREM_U Hodge-number formula (w = 001)
for name, d in [('cyc2_001', 6), ('cyc4_013', 6), ('cyc8_013', 8), ('cyc3_012', 8), ('cyc5_001', 6), ('cyc20_0-1-11', 10)]:
    t0 = time.time()
    R, rng, tries = make_ring(name, d, p)
    P = pieces_of(R)
    tot = sum(x.dim for x in P.values())
    exp = (d - 1) * ((d - 1) ** 3 + 1) // d
    sym = all(P[(a, c)].h[0] == (P[((d - a) % d, (-c) % R.n)].h[2] if ((d - a) % d, (-c) % R.n) in P else 0)
              for (a, c) in P)
    ok = (tot == exp) and sym
    msg = ''
    if name.endswith('_001'):
        n = R.n
        good = all(P[(a, c)].h[0] == c_n(d - 3 - a, c - 1, n) and P[(a, c)].h[2] == c_n(a - 3, -c - 1, n) for (a, c) in P)
        ok = ok and good
        msg = f'THEOREM_U h20/h02 formula: {good}'
    fails += (not ok)
    say(f'[koszul] {name} d{d}: smooth after {tries} redraws (Hilbert = Koszul for e <= 3d-5), sum dim = {tot} vs {exp}, '
        f'h20(a,chi) = h02(-a,-chi): {sym} {msg} -> {"PASS" if ok else "FAIL"} ({time.time()-t0:.1f}s)')

# 2. Taylor coefficients of trivial arcs
for name, d in [('cyc4_013', 6), ('cyc8_013', 8), ('cyc5_012', 10)]:
    R, rng, _ = make_ring(name, d, p)
    P = pieces_of(R)
    K = 4
    arcF = Arc(dict(R.F), K, p)
    dinv = pow(d, p - 2, p)
    okF, okX, cnt = True, True, 0
    for key, x in P.items():
        if x.h[0] == 0:
            continue
        cnt += 1
        W = taylor_W(R, x, arcF, K)
        for k in range(K + 1):
            c = binom_mod((x.a * dinv - 1) % p, k, p)
            if ((W[k] - c * np.vstack([np.eye(x.h[0], dtype=np.int64), np.zeros((x.dim - x.h[0], x.h[0]), np.int64)])) % p).any():
                okF = False
        for i in range(3):
            deltas = []
            for j in range(1, K + 1):
                deltas.append({m: (cf * math.comb(m[i], j)) % p for m, cf in R.F.items() if math.comb(m[i], j)})
            Wg = taylor_general(R, x, deltas, K)
            for col, b in enumerate(x.lev[0]):
                for k in range(K + 1):
                    c = binom_mod((-(b[i] + 1)) % p, k, p)
                    v = np.zeros(x.dim, np.int64); v[col] = c
                    if ((Wg[k][:, col] - v) % p).any():
                        okX = False
    ok = okF and okX
    fails += (not ok)
    say(f'[taylor] {name} d{d}: {cnt} pieces, arc (1+t)F: W_k = binom(a/d-1,k) w_P: {okF}; arcs x_i->(1+t)x_i: '
        f'W_k = binom(-(b_i+1),k) w_P: {okX} -> {"PASS" if ok else "FAIL"}')

# 3. commutant = 1 on THEOREM_U families (non-pure pieces), curve-split pieces > 1 (= 1 + delta^2), hull = dim - delta
for name, d in [('cyc2_001', 6), ('cyc3_001', 7)]:
    R, rng, _ = make_ring(name, d, p)
    P = pieces_of(R)
    vals = [hom_commutant_dim(x, x, p) for x in P.values() if x.h[0] + x.h[2] > 0]
    ok = all(v == 1 for v in vals)
    fails += (not ok)
    say(f'[commutant] {name} d{d}: non-pure pieces {len(vals)}, commutant dims {sorted(set(vals))} -> {"PASS" if ok else "FAIL"}')

for name, d in [('cyc9_013', 7), ('cyc12_015', 8), ('cyc8_016', 9), ('cyc12_014', 9)]:
    R, rng, _ = make_ring(name, d, p)
    n, w = R.n, R.w
    planes = curve_planes(d, n, w, R.S)
    P = pieces_of(R)
    arcs = [Arc({m: rng.randrange(p) for m in R.S}, 4, p) for _ in range(5)]
    lines = []
    ok = True
    for key, x in sorted(P.items()):
        dl = curve_delta(planes, key[0], key[1], d, n)
        if dl == 0:
            continue
        cd = hom_commutant_dim(x, x, p)
        hl = hull_rank([taylor_W(R, x, A, 4) for A in arcs], p) if x.h[0] else None
        nonpure = x.h[0] + x.h[2] > 0
        if nonpure:
            # reducible (E = U + P): first-order commutant must exceed 1; hull must equal dim - delta (h20 > 0)
            good = (cd > 1) and (hl is None or hl == x.dim - dl)
            ok = ok and good
        lines.append(f'{key} h={x.h} delta={dl} comm={cd} hull={hl}')
    # control: non-curve pieces with h20>0 of the same a should have comm 1 (typically)
    fails += (not ok)
    say(f'[curve] {name} d{d}: planes {[(pl["i"], pl["M"], pl["g"]) for pl in planes]}; ' + '; '.join(lines) +
        f' -> {"PASS" if ok else "FAIL"}')

# 4. J-test positive controls + valid perms + gmax
for name, d in [('cyc8_013', 6), ('cyc5_014', 6)]:
    R, rng, _ = make_ring(name, d, p)
    n, w = R.n, R.w
    P = pieces_of(R)
    arcs = [Arc({m: rng.randrange(p) for m in R.S}, 4, p) for _ in range(6)]
    vp = valid_perms(d, n, w, R.S)
    gm = gmax_check(d, n, w, R.S)
    Wc = {}
    def Wof(x):
        if x not in Wc:
            Wc[x] = [taylor_W(R, P[x], A, 4) for A in arcs]
        return Wc[x]
    res = []
    for key, x in sorted(P.items()):
        if x.h[0] >= 1:
            J = jtest(x, x, Wof(key), Wof(key), 4, p)
            res.append((key, J))
            if len(res) >= 4:
                break
    say(f'[jtest] {name} d{d}: J(x,x) = {res}; perms {[(v["pi"], v["ii"], v["iii"]) for v in vp]}; |H|,|Gamma img| = {gm}')
    for v in vp:
        if v['ii'] and v['iii'] and v['order'] == 2:
            imgs = [(k, (k[0], v['cmap'][k])) for k in sorted(P) if v['cmap'].get(k) is not None and v['cmap'][k] != k[1]][:3]
            # J on an automorphism pair
            for (x, y) in imgs:
                if P[x].h[0] >= 1 and y in P:
                    say(f'   automorphism pair {x}->{y}: J = {jtest(P[x], P[y], Wof(x), Wof(y), 4, p)}  (expected >= 1)')
                    break

# 5. L4' on a few h20 = 2 pieces
R, rng, _ = make_ring('cyc3_012', 8, p)
P = pieces_of(R)
out = []
for key, x in sorted(P.items()):
    if x.h[0] == 2:
        out.append((key, x.h, l4p_test(x, random.Random(1)), hom_commutant_dim(x, x, p)))
say('[L4p] cyc3_012 d8 (piece, h, L4p, comm):', out[:8])

say('FAILS', fails)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'outputs', 'unit_tests.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(OUT) + '\n')
