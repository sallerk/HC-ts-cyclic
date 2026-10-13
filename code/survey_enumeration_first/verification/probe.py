import flint, sys
print("flint", flint.__version__, sys.version)
p = 40009
M = flint.nmod_mat(3, 4, [1, 2, 3, 4, 2, 4, 6, 8, 0, 1, 1, 1], p)
R, r = M.rref()
print("rref", R, r, M.rank())
X, nl = M.nullspace()
print("null", X, nl)
print(M.transpose() * M)
print(type(M[0, 1]), int(M[0, 1]))
e = M.entries()
print(type(e), len(e), int(e[5]))
a = flint.nmod_poly([1, 2, 1], p)
b = flint.nmod_poly([1, 1], p)
print(a.gcd(b), a.degree())
Z = flint.nmod_mat(2, 2, p)
print(Z)
try:
    print(M[0:2, 0:2])
except Exception as ex:
    print("noslice", ex)
import numpy
print("numpy", numpy.__version__)
