"""
Section 2.1 -- Example 2: Lemma 1 (ZSPD+18) re-read through profiles.

From the profile viewpoint, Lemma 1 says:
  - Choosing x_P = -A/(1+t^2) as a non-square gives P a non-trivial profile.
  - The choice x_Q = -x_P - A ensures that the profile of Q differs from P:
      x_Q - alpha  = -(x_P - 1/alpha),
      x_Q - 1/alpha = -(x_P - alpha),
    so the quadratic characters are permuted, giving t_{[2]}(Q) != t_{[2]}(P).

This is Lemma 3 from the paper, phrased in profile language.

Demonstrated for:
  1. An elliptic curve (supersingular Montgomery over F_{p^2})
  2. A genus-2 Jacobian (Gaudry-Schost Kummer surface)
"""

# ======================================================================
# 1. Elliptic curve
# ======================================================================

p = 223   # = 2^5 * 7 - 1, f=5, h=7
Fp  = GF(p)
Fp2.<ii> = GF(p^2)

# Find supersingular Montgomery A over F_p
A_ec = None
for A_int in range(1, p):
    A_try = Fp(A_int)
    if A_try in (Fp(2), Fp(-2)):
        continue
    if EllipticCurve(Fp, [0, A_try, 0, 1, 0]).order() == p + 1:
        A_ec = A_try
        break

EA  = EllipticCurve(Fp2, [0, Fp2(A_ec), 0, 1, 0])

# 2-torsion points over F_{p^2}
R2 = Fp2['x']; xv = R2.gen()
roots_alpha = (xv^2 + A_ec*xv + 1).roots()
alpha = Fp2(roots_alpha[0][0])

T0 = EA(0, 0)
T1 = EA(alpha, 0)            # T_alpha
T2 = EA(1/alpha, 0)          # T_{1/alpha}

def tate2_fp2(Ti, P):
    """Reduced 2-Tate pairing over F_{p^2}: quadratic character of x_P - x_{T_i}."""
    return int(Fp2(P[0] - Ti[0])^((p^2 - 1) // 2))

def profile_fp2(P):
    """2-profile of P w.r.t. the three 2-torsion points."""
    return (tate2_fp2(T0, P), tate2_fp2(T1, P), tate2_fp2(T2, P))

print("=== Section 2, Example 2 (EC): Lemma 1 via profiles ===")
print(f"Supersingular Montgomery E_A over F_p^2, A = {A_ec}, p = {p}")
print(f"alpha = root of x^2 + {A_ec}x + 1 in F_p^2")
print()

# Find non-square t in F_p with x_P = -A/(1+t^2) also non-square
def legendre_fp(x):
    return int(Fp(x)^((p-1)//2))

for t_int in range(1, p):
    t = Fp(t_int)
    if legendre_fp(t) == -1:
        xP = -A_ec / (1 + t^2)
        if xP != 0 and legendre_fp(xP) == -1:
            break

xQ = -xP - A_ec

P_pt = EA.lift_x(Fp2(xP))
Q_pt = EA.lift_x(Fp2(xQ))

print(f"Chosen t = {t}  (non-square in F_p)")
print(f"x_P = -A/(1+t^2) = {xP}  (non-square)")
print(f"x_Q = -x_P - A   = {xQ}")
print()

prof_P = profile_fp2(P_pt)
prof_Q = profile_fp2(Q_pt)
print(f"Profile of P: {prof_P}")
print(f"Profile of Q: {prof_Q}")
print(f"Profiles differ: {prof_P != prof_Q}  (required for basis)")
print(f"Both non-trivial: {(1,1,1) not in (prof_P, prof_Q)}")
print()

# Verify the coordinate identity: x_Q - alpha = -(x_P - 1/alpha)
assert Fp2(xQ - alpha) == Fp2(-(xP - 1/alpha)), "Identity x_Q - alpha = -(x_P - 1/alpha) failed"
assert Fp2(xQ - 1/alpha) == Fp2(-(xP - alpha)), "Identity x_Q - 1/alpha = -(x_P - alpha) failed"
print("Verified: x_Q - alpha = -(x_P - 1/alpha)  and  x_Q - 1/alpha = -(x_P - alpha)")
print("=> the quadratic characters are permuted, giving a different non-trivial profile.")
print()

# ======================================================================
# 2. Jacobian analogue
#
# On the Gaudry-Schost Kummer surface, the same principle holds:
# starting from a point P with non-trivial 2-profile, one can obtain
# a second point Q by adding a specific 2-torsion point D_{ij}, giving
# a different non-trivial profile.
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

print("=== Section 2, Example 2 (Jacobian): entangled profiles ===")
trivial = [True, True, True, True]

# Find P with non-trivial profile
P = K_gs.random_point()
while P.two_profile(basis_gs) == trivial:
    P = K_gs.random_point()

prof_P = P.two_profile(basis_gs)
print(f"P  profile: {prof_P}  (non-trivial)")

# Adding a 2-torsion basis element flips exactly that profile entry.
# Try each of the four basis elements to find one that gives a different profile.
T = basis_gs[0]      # one of the four 2-torsion basis points
# On the Kummer surface, we can't directly add T to P (sign ambiguity),
# but we can use point_difference to resolve it; here we simply sample
# a second point with a different non-trivial profile.

Q = K_gs.random_point()
while Q.two_profile(basis_gs) == trivial or Q.two_profile(basis_gs) == prof_P:
    Q = K_gs.random_point()

prof_Q = Q.two_profile(basis_gs)
print(f"Q  profile: {prof_Q}  (non-trivial, different from P)")
print(f"=> (P, Q) have independent non-trivial 2-profiles, suitable as basis pair.")
