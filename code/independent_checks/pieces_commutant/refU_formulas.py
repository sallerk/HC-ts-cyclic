# Referee U: independent check of Sections 2, 3, 6, 7 and (XU).  Own code.
# Usage: python -I -B refU_formulas.py DMAX_HS DMAX_ALL OUTFILE
import sys
from math import gcd

def in_class(d, n):
    return 2 <= n <= d and (d % n == 0 or (d - 1) % n == 0)

def c_brute(n, D, g):
    if D < 0:
        return 0
    return sum(D - b2 + 1 for b2 in range(0, D + 1) if (b2 - g) % n == 0)

def c_closed(n, D, g):
    if D < 0:
        return 0
    r = g % n
    if r > D:
        return 0
    K = (D - r) // n
    return (K + 1) * (D - r + 1) - n * K * (K + 1) // 2

def hilbert(d, n):
    """coefficients of (1-t^{d-1})^2 (1-t^{d-1}s^{-1}) / ((1-t)^2(1-ts)) mod s^n-1"""
    E = 3 * d + 2
    den = [[c_brute(n, e, g) for g in range(n)] for e in range(E)]
    # numerator terms (coef, tdeg, sdeg)
    num = [(1, 0, 0), (-2, d - 1, 0), (1, 2 * d - 2, 0),
           (-1, d - 1, -1), (2, 2 * d - 2, -1), (-1, 3 * d - 3, -1)]
    hs = [[0] * n for _ in range(E)]
    for e in range(E):
        for g in range(n):
            tot = 0
            for (c, td, sd) in num:
                if e - td >= 0:
                    tot += c * den[e - td][(g - sd) % n]
            hs[e][g] = tot
    return hs

def basis_counts(d, n):
    E = 3 * d + 2
    cnt = [[0] * n for _ in range(E)]
    if d % n == 0:
        for i in range(d - 1):
            for j in range(d - 1):
                for k in range(d - 1):
                    cnt[i + j + k][k % n] += 1
    else:
        for j in range(d - 1):
            for i in range(d):
                for k in range(d - 2):
                    cnt[i + j + k][k % n] += 1
            cnt[j + d - 2][(d - 2) % n] += 1
    return cnt

def dimE_formula(d, n, a, chi):
    if d % n == 0:
        return (d - 2) * (d // n - (1 if chi % n == 0 else 0)) + (1 if (chi + a) % n == 0 else 0)
    return (d - 1) ** 2 // n - (d - 2) * (1 if chi % n == 0 else 0)

def units(L):
    return [t for t in range(1, L) if gcd(t, L) == 1] if L > 1 else [1]

def Y_closed(d, n):
    if (d, n) == (6, 2): return 73
    if (d, n) == (6, 3): return 69
    if (d, n) == (6, 6): return 681
    if n == d: return (d - 1) * (d * d - d - 1)
    if 2 * n == d and n >= 4: return n * (d - 1) + d * (d - 2)
    return n * (d - 1)

def C_closed(d, n):
    if (d, n) == (6, 6): return 35
    if n == d - 1: return n - 1
    if n == d: return d - 1
    return n

def dim_list(d, n):
    if (d, n) == (6, 2): return 76
    if (d, n) == (6, 3): return 73
    if (d, n) == (6, 6): return 717
    if n == d - 1: return d * (d - 1)
    if 2 * n == d and d >= 8: return 1 + 3 * d * d // 2 - 2 * d
    if n == d and d >= 7: return 1 + d * (d - 1) ** 2
    generic = (d % n == 0 and 3 * n <= d) or ((d - 1) % n == 0 and 2 * n <= d - 1)
    assert generic, (d, n)
    return 1 + n * d

def pure_list_31(d, n):
    """Prop 3.1 list of pure pieces: dict (a,chi)->dim"""
    out = {}
    if n == d:
        for a in range(1, d):
            for chi in range(1, d):
                if a + chi in (d - 1, d, d + 1):
                    out[(a, chi)] = d - 2 + (1 if a + chi == d else 0)
    if n == d - 1:
        for a in range(1, d):
            out[(a, 0)] = 1
            for chi in range(1, d - 1):
                if a + chi in (d - 1, d):
                    out[(a, chi)] = d - 1
    if 2 * n == d:
        out[(n - 1, 0)] = d - 2; out[(n, 0)] = d - 1; out[(n + 1, 0)] = d - 2
    if 2 * n == d - 1:
        out[(n, 0)] = d; out[(n + 1, 0)] = d
    if (d, n) == (6, 2):
        out[(3, 0)] = 8
    return out

def pure_orbit_list_32(d, n):
    out = set()
    if n == d and d >= 7:
        out = {(a, d - a) for a in range(1, d)}
    elif n == d == 6:
        out = {(a, c) for a in range(1, 6) for c in range(1, 6) if a + c in (5, 6, 7)}
    elif n == d - 1:
        out = {(a, 0) for a in range(1, d)}
    elif 2 * n == d and n >= 4:
        out = {(n, 0)}
    elif (d, n) == (6, 3):
        out = {(2, 0), (3, 0), (4, 0)}
    elif (d, n) == (6, 2):
        out = {(3, 0)}
    return out

def run(d, n, do_hs, do_blocks, log):
    fails = []
    L = d * n // gcd(d, n)
    U = units(L)
    Un = units(n)
    if do_hs:
        hs = hilbert(d, n)
        bc = basis_counts(d, n)
        for e in range(3 * d + 2):
            for g in range(n):
                if hs[e][g] != bc[e][g]:
                    fails.append(("HS!=basis", e, g, hs[e][g], bc[e][g]))
                if 0 <= 3 * d - 6 - e:
                    if hs[e][g] != hs[3 * d - 6 - e][(-2 - g) % n]:
                        fails.append(("gorenstein", e, g))
                elif hs[e][g] != 0:
                    fails.append(("beyond socle", e, g))
        if hs[3 * d - 6][(-2) % n] != 1:
            fails.append(("socle", hs[3 * d - 6]))
    # pieces
    H = {}
    for a in range(1, d):
        for chi in range(n):
            D0 = d - 3 - a
            h20 = c_closed(n, D0, chi - 1)
            h02 = c_closed(n, a - 3, -chi - 1)
            dm = dimE_formula(d, n, a, chi)
            if do_hs:
                lv = [hs[D0 + q * d][(chi - 1) % n] if D0 + q * d >= 0 else 0 for q in range(3)]
                if lv[0] != h20 or lv[2] != h02 or sum(lv) != dm:
                    fails.append(("hodge", a, chi, lv, h20, h02, dm))
            h11 = dm - h20 - h02
            if h11 < 0:
                fails.append(("h11<0", a, chi))
            H[(a, chi)] = (h20, h11, h02)
    # conjugation symmetry
    for (a, chi), h in H.items():
        hc = H[((d - a) % d, (-chi) % n)]
        if hc != (h[2], h[1], h[0]):
            fails.append(("conj", a, chi))
    dim = {s: sum(h) for s, h in H.items()}
    nonzero = [s for s in H if dim[s] > 0]
    zero = sorted(s for s in H if dim[s] == 0)
    exp_zero = sorted((a, 0) for a in range(1, d)) if n == d else []
    if zero != exp_zero:
        fails.append(("zero", zero))
    pure = {s for s in nonzero if H[s][0] == 0 and H[s][2] == 0}
    nonpure = {s for s in nonzero if H[s][0] + H[s][2] > 0}
    pl = pure_list_31(d, n)
    if set(pl) != pure or any(pl[s] != dim[s] for s in pl):
        fails.append(("prop31", sorted(pure), sorted(pl)))
    def orb(s):
        return frozenset(((t * s[0]) % d, (t * s[1]) % n) for t in U)
    orbits = {}
    for s in nonzero:
        orbits.setdefault(orb(s), None)
    pure_orbit_pieces = set()
    for O in orbits:
        for s in O:
            if dim[s] != dim[next(iter(O))]:
                fails.append(("galois dim", s))
        if all(s in pure for s in O):
            pure_orbit_pieces |= set(O)
    if pure_orbit_pieces != pure_orbit_list_32(d, n):
        fails.append(("prop32", sorted(pure_orbit_pieces), sorted(pure_orbit_list_32(d, n))))
    pure_in_np = sorted(pure - pure_orbit_pieces)
    # Separation (a)
    for a in range(1, d):
        seen = {}
        for chi in range(n):
            s = (a, chi)
            if s in nonpure:
                key = (H[s][0], H[s][2])
                if key in seen:
                    fails.append(("separation", s, seen[key]))
                seen[key] = s
    # weight one
    W = {}
    for chi in range(n):
        h10 = c_closed(n, d - 3, chi - 1); h01 = c_closed(n, d - 3, -chi - 1)
        if do_hs:
            if hs[d - 3][(chi - 1) % n] != h10 or hs[2 * d - 3][(chi - 1) % n] != h01:
                fails.append(("w1 hodge", chi))
        W[chi] = (h10, h01)
    wz = [c for c in W if W[c] == (0, 0)]
    exp_wz = [0] if n in (d - 1, d) else []
    if sorted(wz) != exp_wz:
        fails.append(("w1 zero", wz))
    one = sorted(c for c in W if (W[c][0] == 0) != (W[c][1] == 0))
    exp_one = sorted([1, n - 1]) if n == d else []
    if one != exp_one:
        fails.append(("w1 onesided", one))
    two = {c for c in W if W[c][0] > 0 and W[c][1] > 0}
    seen = {}
    for c in two:
        if W[c] in seen:
            fails.append(("w1 separation", c, seen[W[c]]))
        seen[W[c]] = c
    worb = {}
    for c in W:
        if W[c] != (0, 0):
            worb.setdefault(frozenset((t * c) % n for t in Un), None)
    cm = [O for O in worb if all(c not in two for c in O)]
    if (d, n) == (6, 6):
        if cm != [frozenset({1, 5})]:
            fails.append(("cm66", cm))
    elif cm:
        fails.append(("cm", cm))
    # Y by block rules (Prop 6.2), explicit pair loop if requested
    Y = 0
    if do_blocks:
        for a in range(1, d):
            row = [(a, chi) for chi in range(n) if dim[(a, chi)] > 0]
            for s in row:
                for t in row:
                    ps, pt = s in pure_orbit_pieces, t in pure_orbit_pieces
                    if ps and pt:
                        Y += dim[s] * dim[t]
                    elif s == t:
                        Y += 1
                    else:
                        ok = False
                        for tt in U:
                            s2 = ((tt * s[0]) % d, (tt * s[1]) % n)
                            t2 = ((tt * t[0]) % d, (tt * t[1]) % n)
                            np1, np2 = s2 in nonpure, t2 in nonpure
                            if np1 and np2 and (H[s2][0], H[s2][2]) != (H[t2][0], H[t2][2]):
                                ok = True; break
                            if np1 != np2:
                                ok = True; break
                        if not ok:
                            fails.append(("block not closed", s, t))
    else:
        rho = {}
        for s in pure_orbit_pieces:
            rho[s[0]] = rho.get(s[0], 0) + dim[s]
        Y = len([s for s in nonzero if s not in pure_orbit_pieces]) + sum(v * v for v in rho.values())
    C = 0
    for c in W:
        if W[c] == (0, 0):
            continue
        O = frozenset((t * c) % n for t in Un)
        if O in cm:
            continue
        C += 1
    for O in cm:
        for c in O:
            C += (W[c][0] + W[c][1]) ** 2
    if Y != Y_closed(d, n):
        fails.append(("Y", Y, Y_closed(d, n)))
    if C != C_closed(d, n):
        fails.append(("C", C, C_closed(d, n)))
    if 1 + Y + C != dim_list(d, n):
        fails.append(("dim", 1 + Y + C, dim_list(d, n)))
    # (XU)
    sig = {}
    xu_fail = []
    npo = [s for s in nonzero if s not in pure_orbit_pieces]
    for s in npo:
        key = tuple(H[((t * s[0]) % d, (t * s[1]) % n)] for t in U)
        sig.setdefault(key, []).append(s)
    for key, lst in sig.items():
        us = {s[0] for s in lst}
        if len(us) > 1:
            xu_fail.append(sorted(lst))
    log.write("d=%d n=%d L=%d Y=%d C=%d dim=%d pure_orbit=%d pure_in_nonpure=%s XU=%s %s fails=%d\n" % (
        d, n, L, Y, C, 1 + Y + C, len(pure_orbit_pieces), pure_in_np[:6], "holds" if not xu_fail else "FAILS",
        ("ex=" + str(xu_fail[0][:4])) if xu_fail else "", len(fails)))
    for f in fails[:10]:
        log.write("   FAIL %s\n" % (f,))
    return fails, (not xu_fail)

def shioda_identity(d):
    # dim V_a summed over chi with n = d (any n gives the same V_a)
    n = d
    dv = [sum(dimE_formula(d, n, a, chi) for chi in range(n)) for a in range(1, d)]
    lhs = sum(dv[a - 1] * dv[d - a - 1] for a in range(1, d)) + ((d - 1) * (d - 2)) ** 2
    rhs = ((d - 1) ** 6 + (d - 1)) // d
    return lhs == rhs, lhs, rhs

def main():
    dmax_hs = int(sys.argv[1]); dmax = int(sys.argv[2]); out = sys.argv[3]
    with open(out, "w") as log:
        # c_n closed vs brute
        bad = 0
        for n in range(1, 41):
            for D in range(-2, 120):
                for g in range(-3, n + 3):
                    if c_closed(n, D, g) != c_brute(n, D, g):
                        bad += 1
        log.write("c_n closed vs brute (n<=40, D<120): mismatches=%d\n" % bad)
        # monotonicity Lemma 2.4
        bad = 0
        for n in range(1, 41):
            for D in range(0, 120):
                vals = [c_closed(n, D, r) for r in range(0, min(D, n - 1) + 1)]
                if any(vals[i] <= vals[i + 1] for i in range(len(vals) - 1)) or any(v <= 0 for v in vals):
                    bad += 1
                if any(c_closed(n, D, r) != 0 for r in range(D + 1, n)):
                    bad += 1
        log.write("Lemma 2.4 monotonicity: failures=%d\n" % bad)
        for d in range(3, 41):
            ok, l, r = shioda_identity(d)
            if not ok:
                log.write("Shioda dim identity FAILS d=%d %d %d\n" % (d, l, r))
        log.write("Shioda dimension identity checked d=3..40\n")
        tot = 0; totf = 0; xufail = []
        for d in range(6, dmax + 1):
            for n in range(2, d + 1):
                if not in_class(d, n):
                    continue
                f, xu = run(d, n, d <= dmax_hs, d <= dmax_hs, log)
                tot += 1; totf += len(f)
                if not xu:
                    xufail.append((d, n))
        other = [(d, n) for (d, n) in xufail if n != d and 2 * n != d]
        log.write("TOTAL pairs=%d failures=%d\n" % (tot, totf))
        log.write("XU fails at %d pairs; n=d: %d, n=d/2: %d, other: %s\n" % (
            len(xufail), sum(1 for (d, n) in xufail if n == d), sum(1 for (d, n) in xufail if 2 * n == d), other))
        nd_all = all((d, d) in xufail for d in range(6, dmax + 1))
        nh_all = all((d, d // 2) in xufail for d in range(6, dmax + 1, 2))
        log.write("all n=d fail: %s ; all n=d/2 fail: %s\n" % (nd_all, nh_all))

if __name__ == "__main__":
    main()
