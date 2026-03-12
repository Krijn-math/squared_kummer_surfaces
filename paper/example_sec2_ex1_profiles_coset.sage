"""
Section 2.1 -- Example 1: the 2-profile determines the coset of E(F_q)/[2]E(F_q).

The 2-Tate profile map
    t_{[2]} : E(F_q) -> {+1,-1}^m
(where m = rank of E[2](F_q), at most 2 for elliptic curves)
is precisely the coset map for [2]: its image encodes which coset of
E(F_q)/[2]E(F_q) the point P lies in.

Consequence for basis generation: if we want a basis (P,Q) for E[2^f],
we need P and Q in *different* non-trivial cosets, i.e.
    t_{[2]}(P) != t_{[2]}(Q)    and both non-trivial.

Demonstrated for:
  1. An elliptic curve with full rational 2-torsion
  2. A genus-2 Jacobian (Gaudry-Schost Kummer surface)
"""

# ======================================================================
# 1. Elliptic curve
# ======================================================================

p = 101
F = GF(p)
lam1, lam2, lam3 = F(0), F(3), F(7)
E = EllipticCurve(F, [0, -(lam1+lam2+lam3), 0,
                      lam1*lam2 + lam1*lam3 + lam2*lam3, -lam1*lam2*lam3])
T1, T2, T3 = E(lam1, 0), E(lam2, 0), E(lam3, 0)

def tate2(Ti, P, p):
    return int(GF(p)(P[0] - Ti[0])^((p-1)//2))

def two_profile(P):
    if P.is_zero():
        return None
    return (tate2(T1, P, p), tate2(T2, P, p), tate2(T3, P, p))

trivial_ec = (1, 1, 1)

print("=== Section 2, Example 1 (EC): profile determines coset ===")
print(f"E: y^2 = x(x-3)(x-7) over F_{p},  #E = {E.order()}")
print()

# Map every non-identity point to its coset (represented by its profile)
cosets = {}
for pt in E.points():
    if pt.is_zero():
        continue
    prof = two_profile(pt)
    if prof not in cosets:
        cosets[prof] = []
    cosets[prof].append(pt)

print(f"E(F_p)/[2]E(F_p) has {len(cosets)} cosets (should be 4 since E[2] = Z/2 x Z/2):")
for prof, pts in sorted(cosets.items()):
    tag = " <- trivial coset = [2]E(F_p)" if prof == trivial_ec else ""
    print(f"  profile {prof}: {len(pts)} points{tag}")

print()
# Demonstrate: basis (P,Q) for E[2^f] must use points from different cosets
# Find two points with distinct non-trivial profiles
non_trivial = [(prof, pts[0]) for prof, pts in cosets.items() if prof != trivial_ec]
print("Choosing two basis candidates with different non-trivial profiles:")
for prof, pt in non_trivial[:2]:
    print(f"  profile {prof}: representative x = {pt[0]}")

print()
# The isomorphism E(F_p)/[2]E(F_p) -> {+1,-1}^2 ~ Z/2 x Z/2
# (only 2 of the 3 pairing values are independent since their product is fixed)

# ======================================================================
# 2. Jacobian (genus-2): Gaudry-Schost surface
#
# J[2](F_p) = Z_2^4.  The 2-profile is a map
#    t_{[2]} : J(F_p) -> {True,False}^4
# encoding the coset in J(F_p)/[2]J(F_p).
# There are 2^4 = 16 cosets; a basis (P1,P2) for J[2^k] requires P1, P2
# with independent non-trivial profiles.
# ======================================================================

import sys
sys.path.insert(0, '../code')
from kummer_surface import SquaredKummerSurface

p_gs = (1 << 127) - 1
F_gs = GF(p_gs)
_zero_gs = (F_gs(11), F_gs(-22), F_gs(-19), F_gs(-3))
qo_gs = 1809251394333065553571917326471206521441306174399683558571672623546356726339
qt_gs = 1809251394333065553414675955050290598923508843635941313077767297801179626051
K_gs  = SquaredKummerSurface(_zero_gs, jacobian_order=2*qo_gs, twist_order=2*qt_gs)
basis_gs = K_gs.two_torsion_basis()

print("=== Section 2, Example 1 (Jacobian): profile determines coset ===")
print("Gaudry-Schost Kummer over F_{2^127-1}")
print()

trivial_jac = [True, True, True, True]

# Sample several points and show their profiles (= their coset labels)
profiles_seen = {}
print("Sampling random points and recording their 2-profiles (coset labels):")
while len(profiles_seen) < 6:
    P = K_gs.random_point()
    prof = tuple(P.two_profile(basis_gs))
    if prof not in profiles_seen:
        profiles_seen[prof] = True
        tag = " (trivial)" if list(prof) == trivial_jac else ""
        print(f"  profile {prof}{tag}")

print()
print("Two points P1, P2 with independent non-trivial profiles form a valid")
print("basis for J[2^k].  The profile identifies their coset precisely.")

# Verify: 2P always lands in the trivial coset
for _ in range(5):
    P = K_gs.random_point()
    assert (2*P).two_profile(basis_gs) == trivial_jac

print("\nVerified: 2P is in the trivial coset (all 5 random checks passed).")
