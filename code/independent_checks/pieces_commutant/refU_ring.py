# Referee U: check the chain-point Jacobian ring (Lemma 2.2(b)) by exact linear algebra over Z.
# For each degree e: dim R_e = |basis_e|, every monomial minus its claimed normal form lies in J_e,
# and the claimed basis is independent modulo J_e. Also R_e = 0 beyond the socle and the socle monomial.
# Usage: python -I -B refU_ring.py DMIN DMAX OUTFILE
import sys
import flint

def mons(e):
    return [(i, j, e - i - j) for i in range(e + 1) for j in range(e + 1 - i)]

def nf(e, d):
    i, j, k = e
    c = 1
    while True:
        if j >= d - 1 or i >= d:
            return None
        if i >= 1 and k >= d - 2:
            return None
        if i == 0 and k >= d - 1:
            c = -d * c; i = d - 1; k -= d - 1
            continue
        return (c, (i, j, k))

def in_basis(b, d):
    i, j, k = b
    return j <= d - 2 and ((i <= d - 1 and k <= d - 3) or (i == 0 and k == d - 2))

def rank(rows, ncol):
    if not rows:
        return 0
    return flint.fmpz_mat(rows).rank()

def check(d, log):
    gens = [{(d - 1, 0, 0): d, (0, 0, d - 1): 1}, {(0, d - 1, 0): d}, {(1, 0, d - 2): d - 1}]
    ok = True
    for e in range(0, 3 * d - 3):
        M = mons(e); col = {m: c for c, m in enumerate(M)}
        J = []
        if e >= d - 1:
            for mo in mons(e - (d - 1)):
                for g in gens:
                    row = [0] * len(M)
                    for ex, cf in g.items():
                        row[col[(mo[0] + ex[0], mo[1] + ex[1], mo[2] + ex[2])]] += cf
                    J.append(row)
        rJ = rank(J, len(M))
        B = [m for m in M if in_basis(m, d)]
        if len(M) - rJ != len(B):
            log.write("d=%d e=%d dim mismatch %d vs %d\n" % (d, e, len(M) - rJ, len(B))); ok = False
        if B:
            rows = [r[:] for r in J]
            for b in B:
                row = [0] * len(M); row[col[b]] = 1; rows.append(row)
            if rank(rows, len(M)) != rJ + len(B):
                log.write("d=%d e=%d basis dependent\n" % (d, e)); ok = False
        rows = [r[:] for r in J]
        for m in M:
            row = [0] * len(M); row[col[m]] += 1
            r = nf(m, d)
            if r is not None:
                c, b = r
                if not in_basis(b, d):
                    log.write("nf not in basis %s\n" % (b,)); ok = False
                row[col[b]] -= c
            rows.append(row)
        if rank(rows, len(M)) != rJ:
            log.write("d=%d e=%d normal form NOT valid\n" % (d, e)); ok = False
        if e == 3 * d - 6 and B != [(d - 1, d - 2, d - 3)]:
            log.write("socle %s\n" % (B,)); ok = False
        if e > 3 * d - 6 and B:
            ok = False
    log.write("d=%d chain ring check: %s\n" % (d, "OK" if ok else "FAILED"))

def main():
    dmin, dmax, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    with open(out, "w") as log:
        for d in range(dmin, dmax + 1):
            check(d, log)

if __name__ == "__main__":
    main()
