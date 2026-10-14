import sys
from itertools import permutations
sys.stdout.reconfigure(encoding='utf-8')
# Completeness of the survey list under an exact smoothness test (note, Section 9.1 and Section 11, item 6).
# For a set S of monomials of degree d in x0, x1, x2, the general member of span(S) has a smooth zero locus
# if and only if (i) no variable divides every monomial of S, and (ii) for each k, S contains x_k^d or some
# x_k^(d-1) x_j (j != k). [If (i) fails, every member is reducible. Otherwise the base locus is the set of
# coordinate points e_k with x_k^d not in S; Bertini off the base locus; at e_k the linear part of a member is
# spanned by the monomials x_k^(d-1) x_j.]
# We enumerate S_{n,w} for 6 <= d <= 10, 1 <= n <= 2d, 0 <= w1 <= w2 < n, identify sets that differ by a
# permutation of the coordinates, keep those passing (i)-(ii), and compare with the rows of the survey table.
# Usage: python -I -B code/table/smooth_criterion_enum.py ROWS
#   ROWS = the LaTeX rows written by code/table/survey_table2.py ROWS (lines "d & n & w & |S| & m & Y & route ...").

def S_nw(d, n, w1, w2):
    return frozenset((a, b, d - a - b) for a in range(d + 1) for b in range(d + 1 - a)
                     if (w1 * b + w2 * (d - a - b)) % n == 0)

def canon(S):
    return min(tuple(sorted(tuple(m[p[i]] for i in range(3)) for m in S)) for p in permutations(range(3)))

def cond_i(S):
    return not any(all(m[k] > 0 for m in S) for k in range(3))

def cond_ii(S, d):
    return all(any(m[k] in (d, d - 1) for m in S) for k in range(3))

def read_rows(path):
    rows = []
    for line in open(path, encoding='utf-8'):
        if '&' not in line:
            continue
        cells = [c.strip() for c in line.replace(chr(92) * 2, '').split('&')]
        for half in (cells[:7], cells[8:15]):
            if len(half) < 7 or not half[0].isdigit():
                continue
            d, n = int(half[0]), int(half[1])
            w = [int(x) for x in half[2].split(',')]
            rows.append(dict(d=d, n=n, w=w, S=int(half[3]), m=int(half[4]), Y=int(half[5])))
    return rows

rows = read_rows(sys.argv[1])
table = {}
for r in rows:
    S = S_nw(r['d'], r['n'], r['w'][1], r['w'][2])
    if len(S) != r['S']:
        print('|S| mismatch', r)
    key = (r['d'], canon(S))
    if key in table:
        print('duplicate table family', r, table[key])
    table[key] = r
print('table rows: %d, distinct families: %d' % (len(rows), len(table)))

found = {}
excluded = {}
for d in range(6, 11):
    for n in range(1, 2 * d + 1):
        for w1 in range(n):
            for w2 in range(w1, n):
                S = S_nw(d, n, w1, w2)
                if not S:
                    continue
                key = (d, canon(S))
                if key in found or key in excluded:
                    continue
                if not cond_i(S):
                    excluded[key] = 'i'
                elif not cond_ii(S, d):
                    excluded[key] = 'ii'
                else:
                    found[key] = (n, w1, w2, len(S))
for d in range(6, 11):
    print('d=%d: families with smooth general member %d, table %d' % (
        d, sum(1 for k in found if k[0] == d), sum(1 for k in table if k[0] == d)))
print('families with smooth general member: %d' % len(found))
print('excluded (up to permutation): (i) fails %d, (i) holds and (ii) fails %d' % (
    sum(1 for v in excluded.values() if v == 'i'), sum(1 for v in excluded.values() if v == 'ii')))
extra = sorted(k for k in found if k not in table)
missing = [k for k in table if k not in found]
print('smooth but not in the table: %d' % len(extra))
for k in extra:
    print('  d=%d n=%d w=(0,%d,%d) |S|=%d' % ((k[0],) + found[k]))
print('in the table but not smooth by the criterion: %d' % len(missing))
for k in missing:
    print('  ', table[k])
print('moduli m >= 1 for all table rows:', all(r['m'] >= 1 for r in rows))
print('RESULT:', 'complete' if not extra and not missing else 'DIFFERENCE')
