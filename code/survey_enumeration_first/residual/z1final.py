# merge: z1 (V1 rules w/o kernel rules + commutant/hull5/split), ord3 'full' bounds (z1o3_*.json), J-test (jt_*.json)
# into a per-family table and a list of remaining open blocks.  Output: ../outputs/z1_final_table.txt (+ stdout)
import os, sys, json, glob
sys.stdout.reconfigure(encoding='utf-8')
HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = os.path.join(os.path.dirname(HERE), 'outputs')
v1closed = {}
for line in open(os.path.join(OUTD, 'v1rerun_stdout.txt'), encoding='utf-8'):
    if 'FIXED:' in line:
        tag = line.split()[0] + '_' + line.split()[1].replace('=', '')
        v1closed[tag] = 'FIXED: Y=' in line and line.split('FIXED: ')[1].split()[1] == 'closed=True'
lines, openlist, nclosed = [], [], 0
for fn in sorted(glob.glob(os.path.join(OUTD, 'z1o3_cyc*_d*.json')), key=lambda f: (json.load(open(f))['d'], f)):
    D = json.load(open(fn))
    tag = f"{D['name']}_d{D['d']}"
    J = {}
    jf = os.path.join(OUTD, f'jt_{tag}.json')
    if os.path.exists(jf):
        for b in json.load(open(jf))['blocks']:
            J[(b['s'], b['t'])] = b
    lb_t = ub_t = 0
    tools = set()
    for r in D['blocks']:
        lb, ub = r['z'][1], r['z'][2]
        if r.get('open'):
            closers = []
            if 'z_ord3' in r and r['z_ord3'][2] < ub:
                ub = r['z_ord3'][2]
                closers.append('ord3')
            jb = J.get((r['s'], r['t']))
            if jb and jb.get('jbest') is not None and jb['jbest'] < ub:
                ub = jb['jbest']
                closers.append('J')
            if lb == ub:
                tools.update(closers)
            else:
                o = r.get('ord3') or {}
                openlist.append(f"  {tag:18s} {r['s']}->{r['t']} size={r['size']} h={r['hs']}->{r['ht']} kinds={r['ks']},{r['kt']} "
                                f"aut={r['aut']} type={r['open']} lb={lb} ub={ub} (z={r['z'][2]} ord3={o.get('best')} "
                                f"J={jb.get('jbest') if jb else None} Jdims={[t['dims'] for t in jb['tests'][:2]] if jb else None})")
        else:
            rule = r['z'][0]
            if rule == 'split' or rule.endswith('+split'):
                tools.add('split')
        lb_t += r['size'] * lb
        ub_t += r['size'] * ub
    for o in D['orbits']:
        for x, cl in (o.get('member_certs') or {}).items():
            if o['kind'] == 'rigid' and 'hull5' in ' '.join(cl) and not any(('hull3' in ' '.join(o['member_certs'][y]) or 'comm' in o['member_certs'][y]) for y in o['member_certs']):
                tools.add('hull5')
    closed = lb_t == ub_t
    nclosed += closed
    v1c = v1closed.get(tag)
    how = ('V1-rerun' if v1c else ('+'.join(sorted(tools)) or 'z1-rules(comm certs)')) if closed else '-'
    lines.append(f"{tag:18s} Y={lb_t}..{ub_t:<5} {'CLOSED' if closed else 'open  '} by={how:22s} v1rerun_closed={v1c} "
                 f"z1={D['Ylb']}..{D['Yub']} ord3={D['Ylb_ord3']}..{D['Yub_ord3']}")
out = '\n'.join(lines) + f'\n\nclosed: {nclosed} of {len(lines)}\n\nOPEN BLOCKS (after z1 + ord3 + J):\n' + '\n'.join(openlist) + '\n'
open(os.path.join(OUTD, 'z1_final_table.txt'), 'w', encoding='utf-8').write(out)
print(out)
