# Run the Z1 copy of ord3.py (order-k Taylor tests, k <= 3) for a list of (family, d, i) at 2 primes, <= 8 workers.
# Usage: python -I -B ord3drive.py fam:d:i [fam:d:i ...]      (outputs: ../outputs/ord3_*.json, logs ord3log_*.txt)
import os, sys, io, contextlib, traceback, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v1lib import pick_primes, lcm
OUTD = os.path.join(os.path.dirname(HERE), 'outputs')


def one(task):
    name, d, i, p, seed = task
    import ord3
    buf = io.StringIO()
    t0 = time.time()
    try:
        sys.argv = ['ord3.py', name, str(d), str(i), str(p), str(seed), '8', '3']
        with contextlib.redirect_stdout(buf):
            ord3.main()
        status = 'ok'
    except Exception:
        status = 'ERROR ' + traceback.format_exc().replace('\n', ' | ')
    with open(os.path.join(OUTD, f'ord3log_{name}_d{d}_i{i}_p{p}.txt'), 'w', encoding='utf-8') as f:
        f.write(buf.getvalue())
    return f'{name} d={d} i={i} p={p} seed={seed} {status} {time.time() - t0:.1f}s'


if __name__ == '__main__':
    import multiprocessing as mp
    sys.stdout.reconfigure(encoding='utf-8')
    tasks = []
    for arg in sys.argv[1:]:
        name, d, i = arg.split(':')
        d, i = int(d), int(i)
        n = int(name[3:name.index('_')])
        for seed, p in enumerate(pick_primes(lcm(8, d, n), 2, start=50000)):
            tasks.append((name, d, i, p, seed))
    with mp.Pool(min(8, len(tasks)), maxtasksperchild=1) as pool:
        for line in pool.imap_unordered(one, tasks):
            print(line, flush=True)
