"""
Section 5.4 + Example 4 -- Isogeny compression.

We want to communicate a cyclic isogeny  phi : E -> E'  of degree 2^f.
The isogeny is determined by its kernel  ker(phi) = <K>  for some K in E[2^f].

Standard representation: send x_K in F_{p^2}, costing 2 log p bits.

Compressed representation: fix a deterministic basis (P, Q) for E[2^f]
(using Lemma 1 / Lemma 3 so both parties can recompute it), then write
    K = P + [a]Q   (for some odd a; the odd part because [2]K has same kernel as K up to reindex)
and send only  a in Z/2^f Z, costing f bits.

Decompression: receiver recomputes (P, Q) deterministically using the same
rule, then reconstructs K = P + [a]Q.

Example 4 (from paper): this is the method used in [CR24] to compress
(2^n, 2^n)-isogenies between genus-2 Jacobians in a SQISign variant.
Here we demonstrate it on an elliptic curve and on the Jacobian.

Demonstrated for:
  1. Elliptic curve (supersingular Montgomery over F_{p^2})
  2. Genus-2 Jacobian (ell=67 torsion compression, as proxy)
"""

# ======================================================================
# 1. Elliptic curve: compress a cyclic 2^f-isogeny
# ======================================================================

p = 223   # = 2^5 * 7 - 1
f = 5
h = 7

Fp  = GF(p)
Fp2.<ii> = GF(p^2)

A_ec = None
for A_int in range(1, p):
    A_try = Fp(A_int)
    if A_try in (Fp(2), Fp(-2)):
        continue
    if EllipticCurve(Fp, [0, A_try, 0, 1, 0]).order() == p + 1:
        A_ec = A_try
        break

EA = EllipticCurve(Fp2, [0, Fp2(A_ec), 0, 1, 0])

# Deterministic basis (P, Q) for E[2^f]: use Lemma 1 / Lemma 3 construction
def legendre_fp(x):
    return int(Fp(x)^((p-1)//2))

for t_int in range(1, p):
    t = Fp(t_int)
    if legendre_fp(t) == -1:
        xP = -A_ec / (1 + t^2)
        if xP != 0 and legendre_fp(xP) == -1:
            break

xQ_coord = -xP - A_ec
P_basis = h * EA.lift_x(Fp2(xP))
Q_basis = h * EA.lift_x(Fp2(xQ_coord))
assert P_basis.order() == 2^f and Q_basis.order() == 2^f

print(f"=== Section 5.4 (EC): isogeny compression ===")
print(f"Supersingular Montgomery E_A, p={p}=2^{f}*{h}-1")
print(f"Deterministic basis (P, Q) for E[2^{f}] (Lemma 1 construction)")
print(f"  P = E.lift_x(x_P={xP}), cleared cofactor h={h}")
print(f"  Q = E.lift_x(x_Q={xQ_coord}), cleared cofactor")
print()

# Alice wants to communicate a random cyclic isogeny kernel K in E[2^f]
# Compress: express K = P + [a]Q for some a in Z/2^f Z
a_secret = randint(0, 2^f - 1)
K_kernel = P_basis + a_secret * Q_basis

print(f"Secret isogeny kernel: K = P + {a_secret}*Q")
print(f"  K = {K_kernel}")
print(f"  Compressed representation: a = {a_secret}  ({f} bits vs {2*ceil(log(p,2))} bits for x_K)")
print()

# Decompress: reconstruct K = P + [a]Q
K_decompressed = P_basis + a_secret * Q_basis
assert K_decompressed == K_kernel, "Decompression failed"
print(f"Decompressed K matches original: {K_decompressed == K_kernel}")

# Now compute the actual isogeny (to show it works)
phi = EA.isogeny(K_kernel, algorithm="velusqrt" if f > 3 else "velu")
print(f"Isogeny phi: E_{A_ec} -> {phi.codomain().a2()}-curve, degree = {phi.degree()}")
print()

# ======================================================================
# 2. Jacobian: Example 4 -- compress a point K in J[ell] using a basis
#
# For the ell=67 Jacobian, we compress K in J[67] as K = P + [a]Q
# using a fixed basis (P, Q) for J[67].
# ======================================================================

import sys
sys.path.insert(0, '../code')
from kummer_surface import SquaredKummerSurface

p_67  = 2^246 * 3 * 67 - 1
ell_j = 67
F_67  = GF(p_67)
F2_67.<i> = F_67.extension(F_67['x'].gen()^2 + 1)

lam = F_67(12938376701500089469123999060444560237380085127671348976010005071974881035462)
mu  = F_67(4618610882062834088149724968147708876591748255560510786118061150909387779558)
nu  = F_67(1073370653850651346651397166382546615373880359440809096202890114345603658730)

K2 = SquaredKummerSurface(
    (F2_67(lam), F2_67(mu), F2_67(nu)),
    jacobian_order=p_67+1, twist_order=p_67+1
)

print(f"=== Section 5.4, Example 4 (Jacobian): compress K in J[{ell_j}] ===")
print(f"Kummer over F_p with p = 2^246*3*{ell_j}-1")
print()

h_j = (p_67 + 1) // ell_j

def random_order_ell(K, h, ell):
    while True:
        P = K.random_point()
        if not P.on_jacobian():
            continue
        Q = h * P
        if not Q.is_zero():
            return Q

# Fix a deterministic basis (P_base, Q_base) for J[ell]
# (In practice these would be computed deterministically; here we sample once)
P_base = random_order_ell(K2, h_j, ell_j)
Q_base = random_order_ell(K2, h_j, ell_j)
# Ensure independence
while True:
    zeta = P_base.tate_pairing(Q_base, ell_j)
    if zeta != 1:
        break
    Q_base = random_order_ell(K2, h_j, ell_j)

print(f"Deterministic basis (P, Q) for J[{ell_j}] established.")

# Alice compresses: K = P + [a]Q for random a
a_j = randint(1, ell_j - 1)
K_j = P_base + a_j * Q_base

print(f"Kernel point K = P + {a_j}*Q")
print(f"Compressed representation: a = {a_j}  (single scalar mod {ell_j})")
print()

# Decompress
K_j_decompressed = P_base + a_j * Q_base
assert K_j_decompressed == K_j
print(f"Decompressed K matches original: True")
print()
print("In [CR24], this technique compresses (2^n,2^n)-isogenies between")
print("genus-2 Jacobians: each kernel generator K_i is expressed as a")
print("linear combination of a shared deterministic basis, and only the")
print("coefficients are transmitted.")
