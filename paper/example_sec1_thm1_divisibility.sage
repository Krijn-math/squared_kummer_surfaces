"""
Section 1 -- Theorem 1 (Husemoller): divisibility via the 2-Tate pairing.

Theorem 1. Let E : y^2 = (x - lam1)(x - lam2)(x - lam3) over F_q with lam_i in F_q.
Then  P in [2]E  iff  x_P - lam_i  is a square in F_q for i = 1, 2, 3.

Equivalently: the quadratic character of (x_P - lam_i) equals the 2-Tate pairing
t_2(T_i, P)  where T_i = (lam_i, 0) is the i-th rational 2-torsion point.

P in [2]E  iff  the 2-profile  t_{[2]}(P) = (t_2(T_1,P), t_2(T_2,P), t_2(T_3,P))
is trivial (all values = 1).

Demonstrated for:
  1. An elliptic curve E over F_p
  2. A genus-2 Jacobian (Gaudry-Schost Kummer surface over F_{2^127 - 1})
"""

# ======================================================================
# 1. Elliptic curve over F_p
# ======================================================================

p = 101
F = GF(p)

# Curve: y^2 = x(x - 3)(x - 7) = x^3 - 10x^2 + 21x over F_101
# Three rational 2-torsion points: T1=(0,0), T2=(3,0), T3=(7,0)
lam1, lam2, lam3 = F(0), F(3), F(7)
E = EllipticCurve(F, [0, -(lam1+lam2+lam3), 0,
                      lam1*lam2 + lam1*lam3 + lam2*lam3,
                      -lam1*lam2*lam3])
T1, T2, T3 = E(lam1, 0), E(lam2, 0), E(lam3, 0)


def legendre_ec(x, p):
    """Legendre symbol: +1 if x is a nonzero square mod p, -1 if nonsquare."""
    return int(GF(p)(x)^((p - 1) // 2))


def tate2(Ti, P, p):
    """Reduced 2-Tate pairing  t_2(T_i, P) = Legendre(x_P - x_{T_i}, p)."""
    return legendre_ec(P[0] - Ti[0], p)


def two_profile_ec(P):
    """2-profile of P: tuple (t_2(T1,P), t_2(T2,P), t_2(T3,P)) in {+1,-1}^3."""
    return (tate2(T1, P, p), tate2(T2, P, p), tate2(T3, P, p))


print("=== Section 1, Theorem 1 (EC) ===")
print(f"Curve: y^2 = x(x-3)(x-7) over F_{p},  #E = {E.order()}")
print(f"2-torsion points: T1={T1}, T2={T2}, T3={T3}")
print()

# The trivial profile (1,1,1) indicates P in [2]E.
# Verify: for any P, the doubled point 2P always has the trivial profile.
for _ in range(30):
    P = E.random_point()
    if P.is_zero():
        continue
    assert two_profile_ec(2 * P) == (1, 1, 1), "Expected trivial profile for 2P"

print("Verified: t_{[2]}(2P) = (1,1,1) for 30 random non-identity P.")

# Show the four distinct cosets of E(F_p) / [2]E(F_p) via profiles
print("\nFour cosets of E(F_p)/[2]E(F_p) identified by 2-profile:")
seen = {}
for pt in E.points():
    if pt.is_zero():
        continue
    prof = two_profile_ec(pt)
    if prof not in seen:
        seen[prof] = pt
        tag = "(trivial => in [2]E)" if prof == (1, 1, 1) else ""
        print(f"  profile {prof}: e.g. x = {pt[0]}  {tag}")

print()

# ======================================================================
# 2. Genus-2 Jacobian: Gaudry-Schost Kummer surface over F_{2^127 - 1}
#
# J(F_p) ~ Z_2 x Z_2 x Z_2 x Z_2 x Z_r  (r a large prime).
# The 2-Tate pairing profile of P (with respect to the 4-element 2-torsion
# basis) is trivial iff P in [2]J(F_p).
# ======================================================================

import sys
sys.path.insert(0, '../code')
from kummer_surface import SquaredKummerSurface
from kummer_point import SquaredKummerPoint

p_gs = (1 << 127) - 1
F_gs = GF(p_gs)

_zero_gs = (F_gs(11), F_gs(-22), F_gs(-19), F_gs(-3))
qo_gs = 1809251394333065553571917326471206521441306174399683558571672623546356726339
qt_gs = 1809251394333065553414675955050290598923508843635941313077767297801179626051

K_gs = SquaredKummerSurface(_zero_gs, jacobian_order=2*qo_gs, twist_order=2*qt_gs)
basis_gs = K_gs.two_torsion_basis()   # 4-element basis for J[2]

print("=== Section 1, Theorem 1 (Jacobian) ===")
print("Gaudry-Schost Kummer over F_{2^127-1}")
print()

# Trivial profile iff P in [2]J
trivial = [True, True, True, True]

P_gs = K_gs.random_point()
prof_P  = P_gs.two_profile(basis_gs)
prof_2P = (2 * P_gs).two_profile(basis_gs)

print(f"Profile of   P: {prof_P}  (trivial = {prof_P == trivial})")
print(f"Profile of  2P: {prof_2P}  (trivial = {prof_2P == trivial})")

# Verify for several random points
for _ in range(5):
    P_gs = K_gs.random_point()
    assert (2 * P_gs).two_profile(basis_gs) == trivial, "Expected trivial profile for 2P"

print("\nVerified: 2P always has trivial 2-profile on the Jacobian.")
