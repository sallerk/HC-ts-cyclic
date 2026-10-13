# Unit tests for v1lib: Hilbert function (Koszul), trivial arcs (closed-form Taylor coefficients),
# hull / L4' controls on cyc8_017 d8 against the numbers quoted in S2_lemmas.md section 9.
import os, sys, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v1lib import *
from math import comb

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'outputs', 'unit_tests.txt')
log_lines = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    log_lines.append(s)


def pmul(a, b, p):
    out = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            m = madd(m1, m2)
            out[m] = (out.get(m, 0) + c1 * c2) % p
    return {m: c for m, c in out.items() if c}


def pscale(a, s, p):
    return {m: (c * s) % p for m, c in a.items() if (c * s) % p}


def padd(*ps_):
    p = ps_[0]
    out = {}
    for a in ps_[1:]:
        for m, c in a.items():
            out[m] = (out.get(m, 0) + c) % p
    return {m: c for m, c in out.items() if c}


def binom_mod(x, k, p):
    # x in F_p, binom(x, k) = x(x-1)...(x-k+1)/k!
    num, den = 1, 1
    for i in range(k):
        num = num * ((x - i) % p) % p
        den = den * (i + 1) % p
    return num * pow(den, p - 2, p) % p


def family_monomials(d, w, n):
    return [m for m in monos(d) if (w[0] * m[0] + w[1] * m[1] + w[2] * m[2]) % n == 0]


def coord(R, shift, q, c, P):
    d = R.d
    dims = [R.dimR(qq * d - shift, c) for qq in range(1, 4)]
    std = R.std(q * d - shift, c)
    return sum(dims[:q - 1]) + std.index(P)


def arc_tests(name, d, w, n, p, seed):
    rng = random.Random(seed)
    S = family_monomials(d, w, n)
    F = {m: rng.randrange(1, p) for m in S}
    R = Ring(d, w, n, p, F)
    bad = R.check_hilbert()
    log(f'[{name} p={p}] Hilbert/Koszul mismatches (deg 0..3d-5, all weights): {len(bad)}')
    nfail = ntest = 0
    invd = pow(d, p - 2, p)
    for kind in ('S', 'C'):
        shifts = [3 + a for a in range(1, d)] if kind == 'S' else [3]
        qmax = 3 if kind == 'S' else 2
        for shift in shifts:
            a = shift - 3
            for c in range(n):
                dims = [R.dimR(qq * d - shift, c) for qq in range(1, qmax + 1)]
                dimE = sum(dims)
                if dimE == 0:
                    continue
                for q in range(1, qmax + 1):
                    for P in R.std(q * d - shift, c):
                        idx = coord(R, shift, q, c, P)
                        base = {P: 1}
                        # arc (1+t)F
                        x0 = (-q + (a * invd if kind == 'S' else 0)) % p
                        Fk = {(0, 0, 0): 1}
                        for k in range(1, 4):
                            Fk = pmul(Fk, F, p)
                            vec, dd = R.reduce_poly(shift, q + k, c, pmul(base, Fk, p))
                            sc = (-1) ** k * comb(q + k - 1, k)
                            got = [(sc * x) % p for x in (vec + [0] * dimE)[:dimE]]
                            exp = [0] * dimE
                            exp[idx] = binom_mod(x0, k, p)
                            ntest += 1
                            if got != exp:
                                nfail += 1
                        # arcs x_i -> (1+t) x_i
                        for i in range(3):
                            Ls = [None] + [{m: (cf * comb(m[i], k)) % p for m, cf in F.items() if comb(m[i], k) % p and cf} for k in (1, 2, 3)]
                            L1, L2, L3 = Ls[1], Ls[2], Ls[3]

                            def eta(poly, qq):
                                if not poly:
                                    return [0] * dimE
                                v, _ = R.reduce_poly(shift, qq, c, poly)
                                return (v + [0] * dimE)[:dimE]
                            c1 = [(-q * x) % p for x in eta(pmul(base, L1, p), q + 1)]
                            t1 = eta(pmul(base, L2, p), q + 1)
                            t2 = eta(pmul(base, pmul(L1, L1, p), p), q + 2)
                            c2 = [(-q * x + comb(q + 1, 2) * y) % p for x, y in zip(t1, t2)]
                            u1 = eta(pmul(base, L3, p), q + 1)
                            u2 = eta(pmul(base, pmul(L1, L2, p), p), q + 2)
                            u3 = eta(pmul(base, pmul(L1, pmul(L1, L1, p), p), p), q + 3)
                            c3 = [(-q * x + 2 * comb(q + 1, 2) * y - comb(q + 2, 3) * z) % p for x, y, z in zip(u1, u2, u3)]
                            ei = P[i]
                            for k, got in ((1, c1), (2, c2), (3, c3)):
                                exp = [0] * dimE
                                exp[idx] = binom_mod((-(ei + 1)) % p, k, p)
                                ntest += 1
                                if got != exp:
                                    nfail += 1
    log(f'[{name} p={p}] trivial-arc Taylor tests (orders 1-3, all basis classes, surface + curve pieces): {ntest} run, {nfail} failed')
    return R, S


def cyc8_017_controls(p, seed):
    d, w, n = 8, (0, 1, 7), 8
    rng = random.Random(seed)
    for label, S in (('full', family_monomials(d, w, n)),
                     ('Qs', [(8, 0, 0), (4, 2, 2), (0, 4, 4), (0, 8, 0), (0, 0, 8)])):
        F = {m: rng.randrange(1, p) for m in S}
        R = Ring(d, w, n, p, F)
        assert not R.check_hilbert()
        Tk = [{(0, 0, 0)}]
        for k in range(3):
            Tk.append({madd(a, b) for a in Tk[-1] for b in S})
        a = 2
        shift = 3 + a
        for chi in range(8):
            c = chi  # sum of weights = 8 = 0 mod 8
            h = tuple(R.dimR(q * d - shift, c) for q in (1, 2, 3))
            hr, dimE = R.hull(shift, c, S, Tk, 3)
            extra = ''
            if h[0] == 2:
                th = R.theta(shift, c, S, 3)
                M1 = [[th[0][vi][i][0] for vi in range(len(S))] for i in range(h[1])]
                M2 = [[th[0][vi][i][1] for vi in range(len(S))] for i in range(h[1])]
                k = h[1] // 2 + 1
                extra = f' L4cert(k={k})={minrank_cert(M1, M2, k, p, rng)} certified_minrank>={minrank_lower(M1, M2, p, rng)}'
            log(f'[cyc8_017 d8 {label} p={p}] u=2 chi={chi} h={h} hull={hr}{extra}')


if __name__ == '__main__':
    t0 = time.time()
    # primes = 1 mod lcm(8,d,n), not 31081/30241/32003
    pA = pick_primes(lcm(8, 6, 5), 1, start=60000)[0]
    pB = pick_primes(lcm(8, 7, 7), 1, start=60000)[0]
    pC = pick_primes(8, 2, start=60000)
    pD = pick_primes(lcm(8, 6, 2), 1, start=60000)[0]
    arc_tests('cyc5_001 d6', 6, (0, 0, 1), 5, pA, 'T1')
    arc_tests('cyc7_013 d7', 7, (0, 1, 3), 7, pB, 'T2')
    arc_tests('cyc8_017 d8', 8, (0, 1, 7), 8, pC[0], 'T3')
    arc_tests('cyc2_001 d6', 6, (0, 0, 1), 2, pD, 'T4')
    cyc8_017_controls(pC[0], 'T5')
    cyc8_017_controls(pC[1], 'T6')
    log(f'time {time.time() - t0:.1f}s')
    with open(OUT, 'w') as f:
        f.write('\n'.join(log_lines) + '\n')
