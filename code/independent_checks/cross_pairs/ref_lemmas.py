#!/usr/bin/env python3
# Referee checks of Lemma 2.1/2.2 (closed form N*), Lemma 2.5 (floor criterion vs sigma equality), Lemma 2.6,
# and the n = d Hodge numbers against the Fermat/cyclic-cover description of Lemma 3.2 (residue count).
import sys
from math import gcd
from fractions import Fraction as Fr

def cn(D, g, n):
    if D < 0:
        return 0
    r = g % n; tot = 0; b2 = r
    while b2 <= D:
        tot += D - b2 + 1; b2 += n
    return tot

def Nstar(x, n):
    if x <= 0:
        return 0
    e = x % n
    if e == 0:
        e = n
    num = x * x + n * x + e * (n - e)
    assert num % (2 * n) == 0
    return num // (2 * n)

def frac(x):
    return x - (x.numerator // x.denominator)

bad = 0; cnt = 0
for d in range(6, 151):
    for n in range(2, d + 1):
        if d % n and (d - 1) % n:
            continue
        for a in range(1, d):
            for chi in range(n):
                r = chi % n or n
                z = 1 if chi % n == 0 else 0
                sig = a + r
                h20 = cn(d - 3 - a, chi - 1, n); h02 = cn(a - 3, -chi - 1, n)
                cnt += 1
                if h20 != Nstar(d - 1 - sig, n) or h02 != Nstar(sig - 1 - n - z * n, n):
                    bad += 1
print("Lemma 2.1/2.2 closed form: pieces checked", cnt, "mismatches", bad)

# Lemma 2.6
b6 = 0
for N in range(3, 20001):
    if not any(gcd(u, N) == 1 and N / 2 - 2 <= u <= N / 2 for u in range(max(1, N // 2 - 2), N // 2 + 1)):
        b6 += 1
print("Lemma 2.6 N<=20000 failures:", b6)

# Lemma 2.5: floor criterion vs sigma equality, n | d, m >= 2, d <= 48
b5 = 0; c5 = 0
for d in range(6, 49):
    for n in range(2, d // 2 + 1):
        if d % n:
            continue
        m = d // n
        U = [g for g in range(1, d) if gcd(g, d) == 1]
        for a in range(1, d):
            for chi in range(1, n):
                for b in range(1, d):
                    if b == a:
                        continue
                    for chi2 in range(1, n):
                        sig_eq = all((g * a) % d + ((g * chi) % n) == (g * b) % d + ((g * chi2) % n) for g in U)
                        th = Fr(a - b, d); al = Fr(a, d); ga = Fr(chi, n)
                        crit = ((chi2 - chi - (a - b)) % n == 0) and frac(m * th) != 0
                        if crit:
                            for g in U:
                                ph = frac(g * th)
                                if Fr(1, m) <= ph < 1 - Fr(1, m):
                                    crit = False; break
                                if ph < Fr(1, m):
                                    if not (frac(g * al) > ph and frac(g * ga) + m * ph < 1):
                                        crit = False; break
                        c5 += 1
                        if crit != sig_eq:
                            b5 += 1
print("Lemma 2.5 criterion vs sigma equality (d<=48): cases", c5, "mismatches", b5)

# n = d: triple of E_(a,chi) vs H^1(D_G)_k (x) L_(chi,a) by residues:
# h10_k = d-1-k (k in [1,d-1]), h01_k = k-1; L type (1,0) iff chi + a < d.
b3 = 0
for d in range(6, 151):
    for a in range(1, d):
        for chi in range(1, d):
            if (a + chi) % d == 0:
                continue
            h20 = cn(d - 3 - a, chi - 1, d); h02 = cn(a - 3, -chi - 1, d)
            dim = (d - 2) * 1 + (1 if (chi + a) % d == 0 else 0)
            k = (a + chi) % d
            h10, h01 = d - 1 - k, k - 1
            pred = (h10, h01, 0) if a + chi < d else (0, h10, h01)
            if (h20, dim - h20 - h02, h02) != pred:
                b3 += 1
print("Lemma 3.2 Hodge-type consistency (n=d, d<=150) mismatches:", b3)
