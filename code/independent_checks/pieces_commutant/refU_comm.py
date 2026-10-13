# Referee U: independent computation of the full first-order commutant C(E,F0)
# (NO torus/character decomposition: unknowns are all matrix entries of X_q),
# the connectivity graph Gamma(E), generation (Step 1 of Thm 4.9), and the
# chain-point G°-line (Step 3 of Thm 4.5).  Own normal forms.
# Usage: python -I -B refU_comm.py DMIN DMAX OUTFILE [control]
import sys, time, json
from multiprocessing import Pool

P = (1 << 61) - 1

def in_class(d, n):
    return 2 <= n <= d and (d % n == 0 or (d - 1) % n == 0)

def setup(d, n):
    if d % n == 0:
        basis = [(i, j, k) for i in range(d - 1) for j in range(d - 1) for k in range(d - 1)]
        def nf(e):
            i, j, k = e
            if i <= d - 2 and j <= d - 2 and k <= d - 2:
                return (1, e)
            return None
        return 'fermat', basis, nf
    basis = [(i, j, k) for i in range(d) for j in range(d - 1) for k in range(d - 2)] + \
            [(0, j, d - 2) for j in range(d - 1)]
    def nf(e):
        i, j, k = e
        c = 1
        while True:
            if j >= d - 1 or i >= d:
                return None
            if i >= 1 and k >= d - 2:
                return None
            if i == 0 and k >= d - 1:   # x2^{d-1} = -d x0^{d-1}
                c = -d * c; i = d - 1; k -= d - 1
                continue
            return (c, (i, j, k))
    return 'chain', basis, nf

def family(d, n, only_binary=False):
    S = []
    for m2 in range(0, d + 1, n):
        if only_binary and m2 > 0:
            continue
        for m0 in range(0, d - m2 + 1):
            S.append((m0, d - m2 - m0, m2))
    return S

def mult_maps(levels, S, nf):
    """for each q and m: tgt list (coef,t) or None, src list"""
    idx = [{b: r for r, b in enumerate(l)} for l in levels]
    maps = []
    for q in range(len(levels) - 1):
        mq = []
        if levels[q] and levels[q + 1]:
            for m in S:
                tgt = [None] * len(levels[q]); src = [None] * len(levels[q + 1])
                any_ = False
                for c, b in enumerate(levels[q]):
                    r = nf((b[0] + m[0], b[1] + m[1], b[2] + m[2]))
                    if r is None:
                        continue
                    coef, e = r
                    t = idx[q + 1].get(e)
                    if t is None:
                        raise RuntimeError("product left the piece")
                    if src[t] is not None:
                        raise RuntimeError("two sources")
                    tgt[c] = (coef, t); src[t] = (coef, c); any_ = True
                if any_:
                    mq.append((m, tgt, src))
        maps.append(mq)
    return maps

def commutant_dim(levels, maps):
    N = [len(l) for l in levels]
    off = [0]
    for q in range(len(levels)):
        off.append(off[-1] + N[q] * N[q])
    tot = off[-1]
    parent = list(range(tot)); w = [1] * tot; zero = [False] * tot

    def find(x):
        path = []
        while parent[x] != x:
            path.append(x); x = parent[x]
        root = x
        cum = 1
        for y in reversed(path):
            cum = w[y] * cum % P
            w[y] = cum; parent[y] = root
        return root, (w[path[0]] if path else 1)

    for q in range(len(levels) - 1):
        Nq, Nq1 = N[q], N[q + 1]
        for (m, tgt, src) in maps[q]:
            for r in range(Nq1):
                sr = src[r]
                base1 = off[q + 1] + r * Nq1
                for c in range(Nq):
                    tc = tgt[c]
                    if tc is None and sr is None:
                        continue
                    if tc is not None and sr is not None:
                        u = base1 + tc[1]; v = off[q] + sr[1] * Nq + c
                        al = tc[0] % P; be = sr[0] % P
                        ru, wu = find(u); rv, wv = find(v)
                        if ru == rv:
                            if (al * wu - be * wv) % P:
                                zero[ru] = True
                        else:
                            parent[ru] = rv
                            w[ru] = be * wv % P * pow(al * wu % P, P - 2, P) % P
                            zero[rv] = zero[rv] or zero[ru]
                    elif tc is not None:
                        ru, _ = find(base1 + tc[1]); zero[ru] = True
                    else:
                        rv, _ = find(off[q] + sr[1] * Nq + c); zero[rv] = True
    return sum(1 for x in range(tot) if parent[x] == x and not zero[x])

def connected(levels, maps):
    allb = [(q, r) for q in range(len(levels)) for r in range(len(levels[q]))]
    if not allb:
        return True
    id_ = {x: i for i, x in enumerate(allb)}
    par = list(range(len(allb)))
    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for q in range(len(levels) - 1):
        for (m, tgt, src) in maps[q]:
            for c, tc in enumerate(tgt):
                if tc is not None:
                    a, b = f(id_[(q, c)]), f(id_[(q + 1, tc[1])])
                    par[a] = b
    return len({f(i) for i in range(len(allb))}) == 1

def generation_ok(levels, maps):
    for q in range(len(levels) - 1):
        for r in range(len(levels[q + 1])):
            if not any(src[r] is not None and src[r][0] == 1 for (m, tgt, src) in maps[q]):
                return False
    return True

def chain_line_ok(levels, d, D0):
    l0 = levels[0]
    if not l0:
        return None
    k0 = min(b[2] for lv in levels for b in lv if b[2] <= d - 3)
    jstar = D0 - k0
    return sum(1 for b in l0 if b[1] == jstar) == 1

def job(args):
    d, n, control = args
    t0 = time.time()
    kind, basis, nf = setup(d, n)
    # distinct characters check
    if kind == 'fermat':
        chars_ok = True
    else:
        ch = {(b[0] + b[2], b[1], b[2] % (d - 1)) for b in basis}
        chars_ok = len(ch) == len(basis) == (d - 1) ** 3
    bydeg = {}
    for b in basis:
        bydeg.setdefault((sum(b), b[2] % n), []).append(b)
    S = family(d, n)
    if control == "sub2":
        # positive control: the subfamily C[x0,x1,x2^{2n}] (needs 2n | d); pieces split under mu_{2n}
        S = family(d, 2 * n)
    Sb = family(d, n, only_binary=True)
    res = {"d": d, "n": n, "kind": kind, "chars_ok": chars_ok, "w2": [], "w1": []}
    bad = []
    for a in range(1, d):
        D0 = d - 3 - a
        for chi in range(n):
            levels = [sorted(bydeg.get((D0 + q * d, (chi - 1) % n), [])) if D0 + q * d >= 0 else [] for q in range(3)]
            h = tuple(len(l) for l in levels)
            if sum(h) == 0:
                continue
            maps = mult_maps(levels, S, nf)
            dc = commutant_dim(levels, maps)
            conn = connected(levels, maps)
            gen = generation_ok(levels, maps)
            line = chain_line_ok(levels, d, D0) if (kind == 'chain' and h[0] > 0) else None
            nonpure = h[0] + h[2] > 0
            rec = [a, chi, h, dc, conn, gen, line]
            if control and nonpure and h[0] > 0:
                mb = mult_maps(levels, Sb, nf)
                rec.append(commutant_dim(levels, mb))
            res["w2"].append(rec)
            if nonpure and dc != 1:
                bad.append(("w2 comm", a, chi, h, dc))
            if not nonpure and dc != h[1] ** 2:
                bad.append(("w2 pure control", a, chi, h, dc))
            if h[0] > 0 and not (conn and gen):
                bad.append(("w2 conn/gen", a, chi, h, conn, gen))
            if line is False:
                bad.append(("w2 line", a, chi))
    for chi in range(n):
        levels = [sorted(bydeg.get(((q + 1) * d - 3, (chi - 1) % n), [])) for q in range(2)]
        h = tuple(len(l) for l in levels)
        if sum(h) == 0:
            continue
        maps = mult_maps(levels, S, nf)
        dc = commutant_dim(levels, maps)
        conn = connected(levels, maps)
        gen = generation_ok(levels, maps)
        line = chain_line_ok(levels, d, d - 3) if (kind == 'chain' and h[0] > 0) else None
        rec = [chi, h, dc, conn, gen, line]
        res["w1"].append(rec)
        two = h[0] > 0 and h[1] > 0
        if two and (dc != 1 or not conn or not gen or line is False):
            bad.append(("w1", chi, h, dc, conn, gen, line))
        if not two and dc != max(h) ** 2:
            bad.append(("w1 onesided control", chi, h, dc))
    res["bad"] = bad
    res["time"] = round(time.time() - t0, 1)
    return res

def main():
    dmin, dmax, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    control = len(sys.argv) > 4 and sys.argv[4] == "control"
    if len(sys.argv) > 4 and sys.argv[4] == "sub2":
        tasks = [(d, n, "sub2") for d in range(dmax, dmin - 1, -1) for n in range(2, d + 1)
                 if in_class(d, n) and d % (2 * n) == 0]
        with Pool(6) as pool, open(out, "w") as f:
            for res in pool.imap_unordered(job, tasks):
                dims = sorted(r[3] for r in res["w2"] if r[2][0] + r[2][2] > 0)
                f.write(json.dumps({"d": res["d"], "n": res["n"], "sub2_nonpure_dims": dims}) + "\n")
            f.write("DONE\n")
        return
    tasks = [(d, n, control) for d in range(dmax, dmin - 1, -1) for n in range(2, d + 1) if in_class(d, n)]
    tasks.sort(key=lambda t: -(t[0] ** 4) // t[1])
    with Pool(6) as pool, open(out, "w") as f:
        for res in pool.imap_unordered(job, tasks):
            npc = sum(1 for r in res["w2"] if r[2][0] + r[2][2] > 0)
            npc1 = sum(1 for r in res["w2"] if r[2][0] + r[2][2] > 0 and r[3] == 1)
            two = sum(1 for r in res["w1"] if r[1][0] > 0 and r[1][1] > 0)
            two1 = sum(1 for r in res["w1"] if r[1][0] > 0 and r[1][1] > 0 and r[2] == 1)
            # connectivity of Gamma(E) for non-pure pieces with h20 = 0 (informational)
            disc = [(r[0], r[1]) for r in res["w2"] if r[2][0] == 0 and r[2][2] > 0 and not r[4]]
            summ = {"d": res["d"], "n": res["n"], "kind": res["kind"], "chars_ok": res["chars_ok"],
                    "nonpure_w2": npc, "dimC1": npc1, "twosided_w1": two, "dimC1_w1": two1,
                    "h20zero_disconnected": len(disc), "bad": res["bad"][:5], "nbad": len(res["bad"]),
                    "time": res["time"]}
            if control:
                summ["binary_only_dims"] = sorted({r[7] for r in res["w2"] if len(r) > 7})[:10]
            f.write(json.dumps(summ) + "\n"); f.flush()
        f.write("DONE\n")

if __name__ == "__main__":
    main()
