"""
Section 5.3 -- Sampling specific points of order 2^f.

Using Scholten's construction, a supersingular elliptic curve E over F_{p^2}
gives a genus-2 hyperelliptic curve whose Jacobian J satisfies:

    J(F_p) ≅ Z_{(p+1)/2} x Z_{(p+1)/2} x Z_2 x Z_2

We want points with order divisible by 2^f (lying in the Z_{(p+1)/2}^2 part).

The 2-Tate pairing profile classifies points according to which subgroup they lie in:
  - 3 "good" profiles (out of 16):  P has order divisible by 2^f in Z_{(p+1)/2}^2.
  - 9 "ok" profiles:  add one of three 2-torsion points D_{ij} to get a good profile.
  - 4 "bad" profiles:  P lies in the Z_2 x Z_2 part.

Demonstrated for:
  1. Elliptic curve analogue (sampling points with specific 2-adic valuation)
  2. Genus-2 Jacobian (Gaudry-Schost Kummer surface)
"""

# ======================================================================
# 1. Elliptic curve analogue
#
# For a supersingular Montgomery curve E_A over F_{p^2} with p=2^f*h-1,
# sampling x_P as a *non-square* ensures 2^f | ord(P) (Theorem 1 / Lemma 1).
# ======================================================================

p = 223   # = 2^5 * 7 - 1, f=5, h=7
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
n_EA = EA.order()   # = (p+1)^2 = (2^f * h)^2

print(f"=== Section 5.3 (EC): sampling points of order 2^{f} ===")
print(f"Supersingular Montgomery over F_p^2, p={p}=2^{f}*{h}-1")
print(f"#E(F_p^2) = {n_EA} = {factor(n_EA)}")
print()

def legendre_fp(x):
    return int(Fp(x)^((p-1)//2))

# Method: sample x_P in F_p as a non-square.
# This forces P not in [2]E(F_{p^2}) restricted to F_p, ensuring 2^f | ord(P).
good_count = 0
total = 20
for _ in range(total):
    # Sample random x in F_p, check non-square
    while True:
        x_try = Fp.random_element()
        if x_try != 0 and legendre_fp(x_try) == -1:
            break
    try:
        P_try = EA.lift_x(Fp2(x_try))
    except ValueError:
        continue
    # Clear h-cofactor
    P_order_target = h * P_try
    if P_order_target.is_zero():
        continue
    ord_val = valuation(P_order_target.order(), 2)
    if ord_val == f:
        good_count += 1

print(f"Sampled {total} random points with non-square x-coord in F_p:")
print(f"  {good_count}/{total} had order exactly 2^{f} after cofactor clearing.")
print(f"  (A non-square x_P guarantees P not in [2]E, so order is divisible by 2^f.)")
print()

# ======================================================================
# 2. Jacobian: Gaudry-Schost curve
#
# J(F_p) ≅ Z_2 x Z_2 x Z_2 x Z_2 x Z_r  (r = large prime, ~2^126).
# The 16 possible 2-profiles classify points by which part of J(F_p) they lie in.
# For subgroup membership in J[2^k], we need the 2-profile to be in the set
# of "good" profiles.
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

# Canonical "good" profile for J(F_p) ≅ Z_2^4 x Z_r:
# Points in Z_{qo} subgroup have profile = trivial after [2] mapping.
# Actually, the "right" J has order 2*qo and its generator should have non-trivial profile.
# Here: points on the 'Jacobian' branch (on_jacobian()) vs 'Twist' branch.

print("=== Section 5.3 (Jacobian): 2-Tate profiles classify subgroup membership ===")
print("Gaudry-Schost Kummer over F_{2^127-1}")
print()

trivial = [True, True, True, True]
profile_counts = {}
n_samples = 50

for _ in range(n_samples):
    P = K_gs.random_point()
    branch = "Jacobian" if P.on_jacobian() else "Twist"
    prof = tuple(P.two_profile(basis_gs))
    key = (branch, prof)
    profile_counts[key] = profile_counts.get(key, 0) + 1

print(f"Sampled {n_samples} random Kummer points (profiles + branch):")
for (branch, prof), cnt in sorted(profile_counts.items()):
    tag = " (trivial)" if list(prof) == trivial else ""
    print(f"  [{branch}] profile {prof}: {cnt} times{tag}")

print()
print("Interpretation (Section 5.3 of paper):")
print("  The profile tells us which of the 16 cosets of J(F_p)/[2]J(F_p) P lies in.")
print("  For sampling points of order divisible by 2^k, we need a 'good' profile:")
print("  specifically, profiles that correspond to points in the Z_{(p+1)/2}^2 part.")
print()

# The generator 'gen' from the GS example has known order 2*qo.
# Verify its profile is non-trivial (it generates a large-order point).
_gen_raw = (F_gs(1), F_gs(1), F_gs(1), F_gs(78525529738642755703105688163803666634))
gen = 2 * K_gs.point(_gen_raw, validate=False)
prof_gen = gen.two_profile(basis_gs)
print(f"Known Jacobian generator: profile = {prof_gen}  (non-trivial = {prof_gen != trivial})")
