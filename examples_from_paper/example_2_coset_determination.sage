"""
In this example, we showcase the result from [CEMR24]:
that the 2-tate pairing can tell you above which point another point lies.

More precisely, the 2-Tate profile we saw before can only take 4 values
 - (+1, +1, +1) when P is in [2]E(Fp)
 - or one of (+1, -1, -1), (-1, +1, -1) or (-1, -1, +1)
as we need the product of all three individual pairings to be 1 again. 

This precisely matches the four cosets B_i + [2]E(Fp) that E(Fp)/[2]E(Fp)
decomposes into, and each coset contains the points above a certain 2-torsion point.
This all assumes symmetric 2^k-torsion, and hence why it applies on supersingular
Montgomery curves, although a similar thing happens for ordinary curves

Lemma 2. Let E_A : y^2 = x(x-alpha)(x-1/alpha) be a maximal supersingular
Montgomery curve over F_{p^2} with p = 2^f * h - 1.
Let T_0=(0,0), T_1=(alpha,0), T_2=(1/alpha,0) be the three F_{p^2}-rational
2-torsion points.  For P in E[2^f]:

    [2^{f-1}]P = T_i   iff   t_2(T_i, P) = 1  and  t_2(T_j, P) = -1 for j != i.

In other words: the 2-Tate pairings with all three 2-torsion points tell us
exactly which 2-torsion point [2^{f-1}]P equals.

Demonstrated for:
  1. An elliptic curve (supersingular Montgomery curve)
  2. A genus-2 Jacobian (Gaudry-Schost Kummer surface)
"""
from collections import defaultdict

# ======================================================================
# 1. Elliptic curve
# ======================================================================

print()
print("=== Section 1: The 2-Tate pairing determines the coset")
print("=== To show: [2^{f-1}]P = T_i   iff   t_2(T_i, P) = 1  and  t_2(T_j, P) = -1 for j != i")
print()

f = 4
h = 3
p = 2**f * h - 1

Fp  = GF(p)
Fp2.<ii> = GF(p^2)

# Take supersingular Montgomery A over Fp, then raise to Fp2
A = 6

EA = EllipticCurve(Fp2, [0, Fp2(A), 0, 1, 0])
print(f"Supersingular Montgomery E_A over F_p^2, A = {A}, with 2^k-torsion: 2^{f}")

# Compute alpha: root of x^2 + Ax + 1 = 0 in F_{p^2}
R2 = Fp2['x']; xv = R2.gen()
roots_alpha = (xv^2 + A*xv + 1).roots()
alpha = Fp2(roots_alpha[0][0])

T0 = EA(Fp2(0), Fp2(0))
T1 = EA(alpha, 0)
T2 = EA(1/alpha, 0)

def two_profile(P):
    # Tate profile of degree 2, e.g. all tate pairings of P with respect to the 2-torsion on E
    return ( T0.tate_pairing(P, 2, 2), T1.tate_pairing(P, 2, 2), T2.tate_pairing(P, 2, 2) )    

print("\nVerifying Lemma 2: [2^{f-1}]P = T_i  iff  t_2(T_i,P)=1, t_2(T_j,P)=-1 for j!=i")
print()

# We simply take all points of E and divide them into the four cosets, and verify that these have the correct profile

#the naming here is m = minus, p = plus

profiles = defaultdict(list)

for P in EA.points():
    prof = two_profile(P)
    profiles[prof].append(P)

for prof, pts in profiles.items():
    above = Set([ ((p+1) // 2) * P for P in pts])
    print(f"Every point with profile {prof} lies above {above}")
print()


# ======================================================================
# 2. Jacobian analogue (Gaudry-Schost)
#
# J[2] has 4 basis elements.  The 2-Tate profile with respect to the
# 4-element basis encodes which 2-torsion point P lies above
# ======================================================================

print("=== Jacobian (Gaudry-Schost) analogue of Lemma 2 ===")
print("We do the same as above but for the Gaudry Schost Kummer")
print("Eventhough this is not supersingular, the Sylow-2 torsion is symmetric")
print()

import sys
sys.path.insert(0, '../code')
from kummer_surface import SquaredKummerSurface

p = (1 << 127) - 1
F = GF(p)

# we define the Kummer surface by its zero, and feed the orders into the constructor
_zero = (F(11), F(-22), F(-19), F(-3))
qo = 1809251394333065553571917326471206521441306174399683558571672623546356726339
qt = 1809251394333065553414675955050290598923508843635941313077767297801179626051

K = SquaredKummerSurface(_zero, jacobian_order=2*qo, twist_order=2*qt)

# Similar to T1, T2, T3, we take here a basis of the two-torsion
B = K.two_torsion_basis()   # 4-element basis for J[2]

profiles = defaultdict(list)

# I will pick everything from the Jacobian with order 2*qo, for technical reasons
# to do so, I pick a base point Q of this order, and then sample the others to match
while True:
    Q = K.random_point()
    if Q.on_jacobian():
        break

for i in range(1000):
    P = K.random_point(Q)               # this just ensures we are on the same object as Q, e.g. the Jacobian
    prof = tuple(P.two_profile(B))      # the list is otherwise unhashable
    profiles[prof].append(P)

for prof, pts in profiles.items():
    # here, we multiply with qo*qt to make sure we kill all but the last bit of 2-torsion
    above = Set([ qo * P for P in pts])
    print(f"Every point with profile {prof} lies above {above}")
print()

print("Verified: the 2-coset of P is determined by its 2-profile.")
print("Equivalently, when the Sylow-2 torsion is symmetric,"
print("All these points are above the same 2-torsion point!\n")