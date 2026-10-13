import sys, re, os, json
sys.stdout.reconfigure(encoding='utf-8')
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'outputs')
S1 = os.path.join(BASE, 'survey_enumeration_first', 'enumeration', 'summary_v4.txt')
Z1T = os.path.join(BASE, 'survey_second', 'main', 'z3_table.txt')

z = {}
for sub in [os.path.join('survey_second', 'main'), os.path.join('survey_second', 'further_families')]:
    od = os.path.join(BASE, sub)
    for fn in os.listdir(od):
        m = re.match(r'z3_(cyc\S+)_d(\d+)\.json$', fn)
        if m:
            r = json.load(open(os.path.join(od, fn), encoding='utf-8'))
            z[(m.group(1), int(m.group(2)))] = (r['Ylb'], r['Yub'], r['closed'])

# first implementation: S1 closed values, and Z1 values (Z1's Y as recorded in z3_compare / z1 final table)
first = {}
d = None
for line in open(S1, encoding='utf-8'):
    m = re.match(r'== d=(\d+):', line)
    if m:
        d = int(m.group(1)); continue
    if line.startswith('== '):
        d = None; continue
    if d is None:
        continue
    m = re.match(r'(cyc\S+)\s+(D\s+)?(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+.*?\s(YES|no)\s+(YES|no)\s+(\S+)\s+(\S+)', line)
    if m and m.group(10) == 'YES':
        first[(m.group(1), d)] = int(m.group(11))
z1t = os.path.join(BASE, 'survey_enumeration_first', 'residual', 'z1_final_table.txt')
for line in open(z1t, encoding='utf-8'):
    m = re.match(r'\s*(cyc\S+)\s+d=?(\d+)\s.*?Y\s*=?\s*(\d+)(\.\.(\d+))?', line)
    if m and not m.group(4):
        first.setdefault((m.group(1), int(m.group(2))), int(m.group(3)))


def Yformula(d, n):
    # Theorem A: Y(d,n)
    if (d, n) == (6, 2): return 73
    if (d, n) == (6, 3): return 69
    if (d, n) == (6, 6): return 681
    if n == d - 1: return (d - 1) ** 2
    if n == d: return (d - 1) * (d * d - d - 1)
    if 2 * n == d: return n * (d - 1) + d * (d - 2)
    return n * (d - 1)

agree = dis = 0
for k, (lb, ub, cl) in sorted(z.items()):
    if k in first:
        if first[k] == lb and lb == ub:
            agree += 1
        else:
            dis += 1
            print('DIFF', k, first[k], lb, ub)
print('first-impl closed', len(first), 'compared', agree + dis, 'agree', agree, 'disagree', dis)
# Theorem A check
for k, (lb, ub, cl) in sorted(z.items()):
    name, d = k
    m = re.match(r'cyc(\d+)_001$', name)
    if m:
        n = int(m.group(1))
        print('A', k, lb, ub, 'formula', Yformula(d, n), 'OK' if lb == ub == Yformula(d, n) else 'BAD')
