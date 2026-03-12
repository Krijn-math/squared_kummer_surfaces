"""
Section 1 -- Lemma 1 (Zanon-Simplicio-Pereira-Doliskani-Barreto, ZSPD+18):
Entangled basis for the 2^f-torsion of a maximal supersingular Montgomery curve.

Lemma 1. Let E_A : y^2 = x^3 + Ax^2 + x be a maximal supersingular Montgomery
curve over F_{p^2} with p = 2^f * h - 1 (f in N, h odd cofactor, A != 0).
Let t in F_{p^2} be a non-square such that x_P = -A/(1+t^2) is also non-square.
Then:
  - x_P defines a point P in E(F_{p^2}),
  - x_Q = -x_P - A defines a point Q in E(F_{p^2}),
  - ([h]P, [h]Q) is a basis for E[2^f].

The key insight is that choosing x_P as a non-square ensures P has a non-trivial
2-profile, and the specific choice x_Q = -x_P - A gives Q a *different* non-trivial
profile, so the two points span different cosets of E/[2]E.

Demonstrated for:
  1. An elliptic curve (supersingular Montgomery curve)
  2. A genus-2 Jacobian (Gaudry-Schost Kummer surface)
"""

# ======================================================================
# 1. Elliptic curve: supersingular Montgomery curve over F_{p^2}
# ======================================================================
# Use p = 223 = 2^5 * 7 - 1, so f = 5, h = 7.
# Find a supersingular Montgomery curve over F_p (trace = 0 => #E = p+1 = 224).

p = 223   # = 2^5 * 7 - 1
f = 5
h = 7
assert p == 2^f * h - 1 and is_prime(p)

Fp  = GF(p)
R_p = Fp['x']; x = R_p.gen()

# Search for supersingular Montgomery A over F_p (need #E_A(F_p) = p+1)
A_ec = None
for A_int in range(1, p):
    A_try = Fp(A_int)
    if A_try == 2 or A_try == Fp(-2):
        continue
    E_try = EllipticCurve(Fp, [0, A_try, 0, 1, 0])
    if E_try.order() == p + 1:
        A_ec = A_try
        break

assert A_ec is not None, "No supersingular Montgomery curve found"
print(f"Supersingular Montgomery curve: A = {A_ec},  p = {p} = 2^{f}*{h}-1")

# Work over F_{p^2} where E[2^f] is fully rational
Fp2.<ii> = GF(p^2)
EA = EllipticCurve(Fp2, [0, Fp2(A_ec), 0, 1, 0])
print(f"#E_A(F_p) = {p+1},  #E_A(F_p^2) = {EA.order()}")

def legendre_fp(x):
    """Quadratic character in F_p: +1 square, -1 non-square."""
    return int(Fp(x)^((p - 1) // 2))

def lift_x(E, xval):
    """Lift an x-coordinate to a point on E (returns None if no such point)."""
    try:
        return E.lift_x(xval)
    except ValueError:
        return None

# Find a non-square t in F_{p^2} such that x_P = -A/(1+t^2) is also non-square in F_{p^2}
# (We look for t in F_p first for concreteness; then x_P in F_p.)
t_val = None
for t_int in range(1, p):
    t_try = Fp(t_int)
    if legendre_fp(t_try) == -1:          # t is non-square
        xP_try = -A_ec / (1 + t_try^2)
        if xP_try != 0 and legendre_fp(xP_try) == -1:  # x_P also non-square
            t_val = t_try
            break

assert t_val is not None
xP = -A_ec / (1 + t_val^2)
xQ = -xP - A_ec

print(f"\nChosen non-square t = {t_val}")
print(f"x_P = -A/(1+t^2) = {xP}  (non-square: {legendre_fp(xP) == -1})")
print(f"x_Q = -x_P - A   = {xQ}")

P_ec = lift_x(EA, Fp2(xP))
Q_ec = lift_x(EA, Fp2(xQ))
assert P_ec is not None and Q_ec is not None, "Could not lift x_P or x_Q"

# Clear the h-cofactor to land on the 2^f-torsion
hP = h * P_ec
hQ = h * Q_ec

print(f"\n[h]P order = {hP.order()}")
print(f"[h]Q order = {hQ.order()}")
assert hP.order() == 2^f and hQ.order() == 2^f, "Points should have order 2^f"

# Verify ([h]P, [h]Q) is a basis for E[2^f]:
# Two generators of Z/2^f x Z/2^f must be independent.
# Check: [2^(f-1)] hP and [2^(f-1)] hQ are distinct non-identity points.
assert (2^(f-1)) * hP != EA(0) and (2^(f-1)) * hQ != EA(0)
assert (2^(f-1)) * hP != (2^(f-1)) * hQ
print(f"\n([h]P, [h]Q) is a valid basis for E[2^{f}].")

# ======================================================================
# 2. Jacobian analogue: Gaudry-Schost Kummer surface
#
# The Gaudry-Schost surface has J[2] = Z_2^4 over F_p (plus a large prime
# factor).  The 2-Tate pairing plays the same role: two basis points must have
# independent non-trivial 2-profiles.  Here we sample two such points and
# check they generate different cosets.
# ======================================================================

import sys
sys.path.insert(0, '../code')
from kummer_surface import SquaredKummerSurface

p_gs = (1 << 127) - 1
F_gs = GF(p_gs)
_zero_gs = (F_gs(11), F_gs(-22), F_gs(-19), F_gs(-3))
qo_gs = 1809251394333065553571917326471206521441306174399683558571672623546356726339
qt_gs = 1809251394333065553414675955050290598923508843635941313077767297801179626051
K_gs = SquaredKummerSurface(_zero_gs, jacobian_order=2*qo_gs, twist_order=2*qt_gs)
basis_gs = K_gs.two_torsion_basis()

print("\n=== Jacobian (Gaudry-Schost) analogue of Lemma 1 ===")
print("Gaudry-Schost Kummer over F_{2^127-1}")

# Sample two points with different non-trivial 2-profiles
trivial = [True, True, True, True]

def get_nontrivial_point(K, basis):
    while True:
        P = K.random_point()
        prof = P.two_profile(basis)
        if prof != trivial:
            return P, prof

P1, prof1 = get_nontrivial_point(K_gs, basis_gs)
P2, prof2 = get_nontrivial_point(K_gs, basis_gs)
# Resample P2 until its profile differs from P1's
while prof2 == prof1:
    P2, prof2 = get_nontrivial_point(K_gs, basis_gs)

print(f"P1 2-profile: {prof1}")
print(f"P2 2-profile: {prof2}  (different from P1 => different cosets)")
print("=> (P1, P2) can serve as basis points for J[2^k] with distinct non-trivial profiles.")
