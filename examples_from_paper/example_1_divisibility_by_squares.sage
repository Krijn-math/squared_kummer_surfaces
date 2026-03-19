"""
In this example, we showcase a central theorem for divisibility as a Tate pairing result.

Section 1 -- Theorem 1 (Husemoller): divisibility via the 2-Tate pairing.

Theorem 1. Let E : y^2 = (x - lam1)(x - lam2)(x - lam3) over F_q with lam_i in F_q.
Then  P in [2]E  iff  x_P - lam_i  is a square in F_q for i = 1, 2, 3.

Equivalently: the quadratic character of (x_P - lam_i) equals the 2-Tate pairing
t_2(T_i, P)  where T_i = (lam_i, 0) is the i-th rational 2-torsion point.

Hence, rephrased, we can say that P in [2]E  iff  the 2-profile  
t_{[2]}(P) = (t_2(T_1,P), t_2(T_2,P), t_2(T_3,P)) is trivial (all values = 1).

Demonstrated for:
  1. An elliptic curve E over F_p
  2. A genus-2 Jacobian (Gaudry-Schost Kummer surface over F_{2^127 - 1})
"""

# ======================================================================
# 1. Elliptic curve over F_p
# ======================================================================

print()
print("=== Section 1: Husemoller's theorem for elliptic curves")
print("=== To show: Trivial 2-profiles <--> [2]E")
print()

p = 101
F = GF(p)
R.<x> = F[]

lam1 = F.random_element()
lam2 = F.random_element()
lam3 = F.random_element()

# Curve: y^2 = (x - lam1)(x - lam2)(x - lam3)
E = EllipticCurve([0, -(lam1 + lam2 + lam3), 0, lam1*lam2 + lam2*lam3 + lam1*lam3, -lam1*lam2*lam3])

# Three rational 2-torsion points
T1, T2, T3 = E(lam1, 0), E(lam2, 0), E(lam3, 0)

def squares(P):
    # this works as long as P is on the point at infinity or two-torsion
    return ( (P[0] - lam1).is_square(),  (P[0] - lam2).is_square(),  (P[0] - lam3).is_square() )

def divisible(P):
    # returns True of P is divisible by 2, e.g. in [2]E(Fp)
    if len(P.division_points(2)) > 0:
        return True
    
    return False

# Husemoller's theorem essentially says P has no division point if squares(P) returns [false, false, false]
for i in range(100):
    P = E.random_element()
    
    # we skip the exceptional points for which squares fails
    if P == E(0) or P == T1 or P == T2 or P == T3:
        continue
    
    if divisible(P):
        assert squares(P) == (True, True, True)
    else:
        assert squares(P) != (True, True, True)

print("Husemoller's theorem verified")

# The key idea is that we should interpret Husemoller's theorem as a statement about the Tate pairing
# and in particular the profile. The squares-values are precisely the 2-Tate pairings

def profile(P):
    # Tate profile of degree 2, e.g. all tate pairings of P with respect to the 2-torsion on E
    return ( T1.tate_pairing(P, 2, 1), T2.tate_pairing(P, 2, 1), T3.tate_pairing(P, 2, 1) )

for i in range(100):
    # we skip the exceptional points for which squares fails
    if P == E(0) or P == T1 or P == T2 or P == T3:
        continue
    
    # we should translate True <--> 1 and False <--> -1, as interpretation of the Tate pairing and .is_square
    assert (profile(P) == (1,1,1)) == (squares(P) == (True, True, True) )
    
print("Husemoller's theorem are 2-Tate profiles!")


# and so we can finalise our rephrasing:
# The 2-tate profile of P ∈ E(Fp) is trivial if and only if P ∈ [2]E(Fp)

trivial_profiles = { P for P in E.points() if profile(P) == (1,1,1) }
two_E = {2*P for P in E.points()}

print(f"trivial profiles equal [2]E: {trivial_profiles == two_E}\n")


# ======================================================================
# Now we give a second example of this approch, but over dimension 2 Jacobians
# 
# We use the example of the Gaudry-Schost Kummer surface over F_{2^127 - 1}
# which has torsion subgroup
#      J(F_p) ~ Z_2 x Z_2 x Z_2 x Z_2 x Z_r  (r a large prime).
# The 2-Tate pairing profile of P (with respect to the 4-element 2-torsion
# basis) is trivial iff P in [2]J(F_p).
# ======================================================================


print()
print("=== Section 1: Husemoller's theorem for Jacobians")
print("=== using Gaudry-Schost Kummer over F_{2^127-1}")
print()

import sys
sys.path.insert(0, '../code')
from kummer_surface import SquaredKummerSurface
from kummer_point import SquaredKummerPoint

p = (1 << 127) - 1
F = GF(p)

# we define the Kummer surface by its zero, and feed the orders into the constructor
_zero = (F(11), F(-22), F(-19), F(-3))
qo = 1809251394333065553571917326471206521441306174399683558571672623546356726339
qt = 1809251394333065553414675955050290598923508843635941313077767297801179626051

K = SquaredKummerSurface(_zero, jacobian_order=2*qo, twist_order=2*qt)

# Similar to T1, T2, T3, we take here a basis of the two-torsion
B = K.two_torsion_basis()   # 4-element basis for J[2]

# Trivial profile iff P in [2]J
trivial = [True, True, True, True]

# now if we take a random point P, it might either be in
# J \ [2]J or in [2]J
# in any case, [2]P must be in [2]J
for i in range(100):
    P = K.random_point()
    assert (2*P).two_profile(B) == trivial

    # if P itself has trivial profile, we verify that indeed r*P is zero, which implies indeed P in [2]J
    # but r is either qo or qt, so lets do both
    if P.two_profile(B) == trivial:
        assert (qo*qt*P).is_zero()

print("Verified: P has trivial 2-profile if and only if P in [2]J\n")
