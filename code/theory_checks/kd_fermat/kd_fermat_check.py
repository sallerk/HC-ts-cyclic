import sys
from math import gcd
from collections import Counter
sys.stdout.reconfigure(encoding='utf-8')
# For n = d >= 7: k(d) = |KR(d)| (Definition 7.5) against rho(X_d) - 2d^2 + 5d - 4, where X_d is the Fermat
# surface of degree d and rho(X_d) = 1 + #B_d (Shioda), B_d = {(a0,a1,a2,a3): a_i != 0 mod d, sum = 0,
# sum_i (g a_i)_d = 2d for all units g}. Also against the values listed in Section 11(5) of the note, and
# Shioda's criterion B_d = D_d (decomposable characters) iff gcd(d,6) = 1 for d >= 7 (note, Remark 10.3).
# KR(d) is counted directly from the definition. The identity is claimed for d >= 7 only; for d = 6 the line
# prints both numbers (k(6) = 16, rho(X_6) - 2*36 + 30 - 4 = 40). Usage: python -I -B kd_fermat_check.py > OUT
LISTED = {6: 16, 8: 84, 9: 96, 10: 208, 12: 412, 14: 480, 15: 456, 16: 436, 18: 1096, 20: 1284, 21: 792,
          22: 880, 24: 2044, 26: 1152, 27: 816, 28: 1540, 30: 3976, 32: 1524, 33: 1200, 34: 1792, 36: 3100,
          38: 2160, 39: 1656, 40: 3220, 42: 7096, 44: 2724, 45: 2376, 46: 2992, 48: 4828, 50: 3408, 51: 2784,
          52: 3652, 54: 4696, 56: 4404, 57: 3456, 58: 4480, 60: 13132, 62: 5040, 63: 4488, 64: 5236, 66: 6760,
          68: 5892, 69: 5016, 70: 6256, 72: 9148, 74: 6912, 75: 6096, 76: 7204, 78: 8872, 80: 8820, 81: 6864,
          82: 8320, 84: 15484, 86: 9072, 87: 7896, 88: 9364, 90: 13864, 92: 10212, 93: 9000, 94: 10672,
          96: 13084, 98: 11568, 99: 10176, 100: 12484}

def rep(x, d):
    r = x % d
    return r if r else d

def kr_direct(d):
    units = [g for g in range(1, d) if gcd(g, d) == 1]
    sig = {}
    for a in range(1, d):
        for c in range(1, d):
            if a + c == d or (d == 6 and a + c in (5, 7)):
                continue  # pieces in pure orbits: sigma = d for d >= 7, sigma in {5,6,7} for d = 6
            sig[(a, c)] = tuple(rep(g * a, d) + rep(g * c, d) for g in units)
    groups = {}
    for s, v in sig.items():
        groups.setdefault(v, []).append(s)
    k = 0
    for v, mem in groups.items():
        for s in mem:
            for t in mem:
                if s[0] != t[0]:
                    k += 1
    return k

def hodge_chars(d):
    units = [g for g in range(1, d) if gcd(g, d) == 1]
    cnt = Counter()
    vecs = {}
    for x in range(1, d):
        for y in range(1, d):
            v = tuple(rep(g * x, d) + rep(g * y, d) for g in units)
            vecs[(x, y)] = v
            cnt[v] += 1
    total = 0
    for p, v in vecs.items():
        total += cnt[tuple(2 * d - e for e in v)]
    return total

bad = 0
for d in range(6, 101):
    B = hodge_chars(d)
    rho = 1 + B
    formula = rho - 2 * d * d + 5 * d - 4
    k = kr_direct(d)
    ref = LISTED.get(d, (d - 1) * (d - 3) if gcd(d, 6) == 1 else None)
    lines = 3 * (d - 1) * (d - 2) + (1 if d % 2 == 0 else 0)
    ok = (k == ref) and (d < 7 or k == formula)
    if not ok:
        bad += 1
    if d <= 12 or not ok:
        print('d=%d rho=%d #D=%d k_direct=%d k_listed=%s rho-2d^2+5d-4=%d %s' % (d, rho, lines, k, ref, formula, 'ok' if ok else 'MISMATCH'))
    # Shioda: B_d = D_d iff d <= 4 or gcd(d,6) = 1; with the identity, KR(d) = swap pairs iff B_d = D_d
    if d >= 7 and ((B == lines) != (gcd(d, 6) == 1)):
        print('d=%d: B=D test disagrees with gcd(d,6)=1' % d)
        bad += 1
print('checked 6 <= d <= 100; mismatches (d >= 7, or listed values):', bad)
