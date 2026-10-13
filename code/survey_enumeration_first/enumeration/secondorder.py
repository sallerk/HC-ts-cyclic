"""Orders 1-3 (D1's ord3.py, local copy) on selected R3-residual blocks, 2 primes (the survey's primes for the family).
Block in ord3 notation: sigma = monomial G-character of the source H^{2,0} (= survey chi - w.(1,1,1) mod n), Delta = chi' - chi.
Usage: python -I -B secondorder.py TAG SPEC [SPEC ...]   SPEC = family:d:i:sigma:Delta[,sigma:Delta...]"""
import json
import os
import subprocess
import sys

sys.stdout.reconfigure(encoding='utf-8')
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, '..', 'outputs'))
tag = sys.argv[1]
for spec in sys.argv[2:]:
    fam, d, i, blks = spec.split(':', 3)
    d, i = int(d), int(i)
    blks = [tuple(int(x) for x in b.split(':')) for b in blks.split(',')]
    surv = json.load(open(os.path.join(OUT, 'survey_d%d%s.json' % (d, tag)), encoding='utf-8'))
    rec = [r for r in surv if r['name'] == fam]
    primes = rec[0]['primes'] if rec else [30097, 30169]
    for seed, p in enumerate(primes):
        cmd = [sys.executable, '-I', '-B', os.path.join(HERE, 'ord3.py'), fam, str(d), str(i), str(p), str(seed), '8', '3']
        out = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8').stdout
        fn = os.path.join(OUT, 'ord3_%s_d%d_i%d_p%d_s%d_a8_o3.json' % (fam, d, i, p, seed))
        res = json.load(open(fn, encoding='utf-8'))
        print("%s d=%d i=%d p=%d: %s" % (fam, d, i, p, out.splitlines()[0][:160] if out else 'NO OUTPUT'))
        for (sg, De) in blks:
            b = [x for x in res['blocks'] if x['sigma'] == sg and x['Delta'] == De]
            if not b:
                print("   block sigma=%d Delta=%d: not present (no H^{2,0} on one side)" % (sg, De))
                continue
            b = b[0]
            labs = [k for k in ('o1', 'V1o2', 'V1o3', 'V3o2', 'V3o3') if k in b]
            print("   block sigma=%d Delta=%d c%s phibar%s | c-proj %s | full %s" % (
                sg, De, b['c_shape'], b['phibar_shape'], " ".join("%s=%d" % (k, b[k]['c']) for k in labs),
                " ".join("%s=%d" % (k, b[k]['full']) for k in labs)))
