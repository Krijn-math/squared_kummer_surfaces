"""
Section 5.7 -- Subgroup membership testing via Tate pairing profiles.

Let E/F_p be an elliptic curve with a large prime-order subgroup G = E[r](F_p),
and suppose E(F_p) ≅ Z_e^2 x Z_r  for some cofactor e (often e is a small power of 2).

Lemma 6. P in G  iff  the e-Tate profile  t_{[e]}(P)  is trivial.

Proof sketch: P in G = E[r](F_p) = [e]E(F_p)  (since G has prime order r and
[e]P = O iff r | e*ord(P), so [e]E(F_p) = G when gcd(e,r)=1 and r is prime).
Triviality of t_{[e]}(P) is exactly the condition P in [e]E(F_p).

Application: to test P in G, compute the 2-Tate profile (or e-Tate profile)
and check if it is trivial.  This can be much cheaper than computing [r]P.

Remark 3 (FourQ curve): #E(F_{p^2}) = 2^3 * 7^2 * r, G = E[r](F_{p^2}).
Testing P in G via profile t_{[7]}(P) is possible using rational E[7].

Gaudry-Schost remark (Section 5.7): #J(F_p) ≅ Z_2^4 * Z_r.
Testing P in J[r](F_p) via 2-Tate profile with respect to J[2] basis.

Demonstrated for:
  1. Elliptic curve with small cofactor (2-torsion test)
  2. Genus-2 Jacobian (Gaudry-Schost, test P in J[r])
"""

# ======================================================================
# 1. Elliptic curve: membership in prime-order subgroup
# ======================================================================

# Find an elliptic curve over F_p with group structure Z_e^2 x Z_r,
# where e is small and r is a large prime.
# Simple case: e=1 (prime order curve), then E(F_p) = Z_r, and testing
# P in G is trivial.  More interesting: e=2, so E(F_p) = Z_4 x Z_r or Z_2^2 x Z_r.

# Use p=101, E: y^2 = x(x-3)(x-7) (full 2-torsion), #E = some n = 4*k or 2^m*r
p_ec = 101
F_ec = GF(p_ec)
lam1, lam2, lam3 = F_ec(0), F_ec(3), F_ec(7)
E_ec = EllipticCurve(F_ec, [0, -(lam1+lam2+lam3), 0,
                             lam1*lam2+lam1*lam3+lam2*lam3, -lam1*lam2*lam3])
n_ec = E_ec.order()
print(f"=== Section 5.7 (EC): subgroup membership testing ===")
print(f"E: y^2 = x(x-3)(x-7) over F_{p_ec},  #E = {n_ec} = {factor(n_ec)}")

# Factor n and identify the large prime part
fac_n = list(factor(n_ec))
r_ec = max(q^e for q,e in fac_n if is_prime(q))   # largest prime power
e_ec = n_ec // r_ec
print(f"  Large prime component r = {r_ec},  cofactor e = {e_ec}")

T1_ec = E_ec(lam1, 0)
T2_ec = E_ec(lam2, 0)

def tate2(Ti, P):
    return int(F_ec(P[0] - Ti[0])^((p_ec-1)//2))

def two_profile_ec(P):
    return (tate2(T1_ec, P), tate2(T2_ec, P))

trivial_ec = (1, 1)

# The prime-order subgroup G = [e]E(F_p) has order r.
# Membership test: P in G  iff  t_{[e]}(P) is trivial
# (which here means t_2(T1, P) = t_2(T2, P) = 1 when e = 4, or adjust by cofactor)

# Generate a random P in G and check its profile
G_points = [e_ec * pt for pt in E_ec.points() if not pt.is_zero()]
G_points = [pt for pt in G_points if not pt.is_zero()]

print(f"\n  Prime-order subgroup G = [e]E(F_p) has {len(set(G_points))} distinct elements (order r={r_ec}).")
print(f"  For P in G: t_{{[2]}}(P) should be trivial (all +1).")

in_G_wrong = 0
total_in_G = 0
for _ in range(20):
    P = E_ec.random_point()
    in_G = (r_ec * P == E_ec(0))   # direct membership test
    prof = two_profile_ec(P)
    trivial_prof = (prof == trivial_ec)
    if in_G != trivial_prof:
        in_G_wrong += 1
    total_in_G += (1 if in_G else 0)

print(f"\n  20 random points tested: {total_in_G} in G.")
if in_G_wrong == 0:
    print(f"  All {20-in_G_wrong}/20 agreed: P in G  iff  t_{{[2]}}(P) = (1,1).  ✓")
else:
    print(f"  Warning: {in_G_wrong} disagreements (may need e-profile, not 2-profile).")

print()

# ======================================================================
# 2. Jacobian: Gaudry-Schost membership test
#
# J(F_p) ≅ Z_2 x Z_2 x Z_2 x Z_2 x Z_r  (r = 126-bit prime).
# Prime-order subgroup G = J[r](F_p).
# P in G  iff  t_{[2]}(P) is trivial (all four 2-Tate pairings = 1).
# ======================================================================

import sys
sys.path.insert(0, '../code')
from kummer_surface import SquaredKummerSurface

p_gs = (1 << 127) - 1
F_gs = GF(p_gs)
_zero_gs = (F_gs(11), F_gs(-22), F_gs(-19), F_gs(-3))
qo_gs = 1809251394333065553571917326471206521441306174399683558571672623546356726339
qt_gs = 1809251394333065553414675955050290598923508843635941313077767297801179626051

# jacobian_order = 2*qo_gs, twist_order = 2*qt_gs
# The Jacobian has order 2*qo_gs where qo_gs is a large prime.
K_gs  = SquaredKummerSurface(_zero_gs, jacobian_order=2*qo_gs, twist_order=2*qt_gs)
basis_gs = K_gs.two_torsion_basis()

print("=== Section 5.7 (Jacobian): membership in J[r](F_p) via 2-profile ===")
print("Gaudry-Schost Kummer: J(F_p) ≅ Z_2 x Z_2 x Z_2 x Z_2 x Z_r")
print(f"Prime order r = qo_gs (a {len(str(qo_gs))}-digit prime)")
print()

trivial = [True, True, True, True]

# Points in J[r](F_p) are those with [2]*P = O (since order divides 2*qo, and r | qo)
# Actually: on the Kummer surface with order 2*qo, points on the Jacobian have order dividing 2*qo.
# A point is in the prime-order part iff [2]P has order qo_gs... let me be precise:
# J(F_p) ≅ Z_2 * Z_r, so 2-Tate profile trivial <=> P in [2]J(F_p) = J[r](F_p).

print("Testing: P in J[r](F_p)  iff  2-profile t_{[2]}(P) is trivial.")
print()

# Sample points and classify
results = []
for _ in range(10):
    P = K_gs.random_point()
    is_on_jac = P.on_jacobian()
    prof = P.two_profile(basis_gs)
    trivial_prof = (prof == trivial)
    # Direct test: P in J[r](F_p) iff [2]P has order qo (or we check [qo_gs]P == 0)
    in_prime_subgroup = (qo_gs * P).is_zero()
    results.append((trivial_prof, in_prime_subgroup, is_on_jac))
    print(f"  P on_jacobian={is_on_jac},  trivial profile={trivial_prof},  "
          f"P in J[r](F_p)={in_prime_subgroup},  agree={trivial_prof == in_prime_subgroup}")

n_agree = sum(1 for t, d, _ in results if t == d)
print(f"\n{n_agree}/10 profile tests agreed with direct membership test.")
print()
print("This is the key result of Koshelev [Kos23] and Sec. 5.7:")
print("  Subgroup membership testing = checking triviality of Tate profile.")
print("  Cost: 4 cubical pairings of degree 2 (much cheaper than scalar mult by r).")
