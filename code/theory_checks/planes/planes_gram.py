import sys
sys.stdout.reconfigure(encoding='utf-8')
# Gram matrix of h^2 and the n*d planes P_{eta,gamma} = {y = eta*gamma(x)} on X_F = {F(x)+F(y)=0} in P^5,
# F in C[x0,x1,x2^n]_d. Planes indexed by (i, j): eta = exp(pi*i*(2i+1)/d), gamma: x2 -> exp(2 pi i j/n) x2.
# Intersection numbers (clean intersections of linear spaces, eps_j = (1-(1-d)^j)/d):
#   same plane: eps_3; same eta, different gamma: line, eps_2; different eta with eta*e = eta'*e': point, eps_1;
#   otherwise empty: 0. h^2.P = 1, h^4 = d. Rank mod a large prime.
P = 2**61 - 1

def eps(j, d):
    return (1 - (1 - d) ** j) // d

def rank_mod(M):
    M = [[x % P for x in row] for row in M]
    r = 0
    rows, cols = len(M), len(M[0])
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i][c]:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], P - 2, P)
        M[r] = [(x * inv) % P for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % P for a, b in zip(M[i], M[r])]
        r += 1
    return r

def dimformula(d, n):
    if (d, n) == (6, 2): return 76
    if (d, n) == (6, 3): return 73
    if (d, n) == (6, 6): return 717
    if n == d - 1: return d * (d - 1)
    if n == d: return 1 + d * (d - 1) ** 2
    if 2 * n == d: return 1 + 3 * d * d // 2 - 2 * d
    return 1 + n * d

# exact arithmetic on angles: eta = exp(2 pi i (2i+1)/(2d)), gamma scalar exp(2 pi i j / n);
# eta*e as a fraction mod 1: (2i+1)/(2d) + j/n
from fractions import Fraction
bad = 0
for d in range(6, 21):
    for n in range(2, d + 1):
        if d % n and (d - 1) % n:
            continue
        planes = [(i, j) for i in range(d) for j in range(n)]
        def ang(p):
            return (Fraction(2 * p[0] + 1, 2 * d) + Fraction(p[1], n)) % 1
        N = len(planes) + 1
        G = [[0] * N for _ in range(N)]
        G[0][0] = d
        for a in range(1, N):
            G[0][a] = G[a][0] = 1
        for a, p in enumerate(planes, 1):
            for b, q in enumerate(planes, 1):
                if p == q:
                    v = eps(3, d)
                elif p[0] == q[0]:
                    v = eps(2, d)
                elif ang(p) == ang(q):
                    v = eps(1, d)
                else:
                    v = 0
                G[a][b] = v
        r = rank_mod(G)
        rp = rank_mod([row[1:] for row in G[1:]])
        df = dimformula(d, n)
        if n == d:
            exp_r = d * d - d + 1
        elif n == d - 1:
            exp_r = n * d
        else:
            exp_r = n * d + 1
        ok = (r == exp_r)
        generic = (n == d - 1) or ((d - 1) % n == 0) or (d % n == 0 and 3 * n <= d and (d, n) != (6, 2))
        spans = (r == df)
        if not ok or (generic and not spans) or (not generic and spans):
            bad += 1
        print(d, n, 'rank', r, 'planes-only', rp, 'expected', exp_r, 'dimHdg', df, 'spans' if spans else '', '' if ok else 'MISMATCH')
print('bad', bad)
