"""
Section 5.5 -- Computing the ell-Sylow torsion using Tate pairing profiles.

Algorithm (from Section 5.5 of the paper):

Input:  A basis B = {P_1,...,P_r} of A[ell](F_q),
        a division algorithm  Div(P) -> Q  such that [ell]Q = P (if it exists).

Output: A basis B' = {Q_1,...,Q_r} of the Sylow-ell subgroup S_{ell,q}(A).

Steps:
  1. Compute Q_i = P_i  for each i (starting points equal to A[ell] basis).
  2. For each i: apply the division algorithm repeatedly:
       Q_i <- (1/ell) * Q_i   until the profile t_{[ell]}(Q_i) becomes non-trivial.
  3. If profiles of {Q_1,...,Q_r} span mu_ell^r, output B' = {Q_1,...,Q_r}.
     Otherwise, fix linear dependencies:
       - Find j,k with  t_{[ell]}(Q_j) = prod_i t_{[ell]}(Q_i)^{lambda_i}
       - Set Q_j <- Q_j - sum lambda_i * Q_i  and repeat step 2.

The key insight: t_{[ell]}(Q_i) is trivial iff Q_i in [ell]A(F_q), which
means Q_i can be divided further.  The algorithm terminates because the
Sylow subgroup is finite.

Demonstrated for:
  1. Elliptic curve (small example with explicit 2-Sylow computation)
  2. Genus-2 Jacobian (ell=67, Kummer surface)
"""

# ======================================================================
# Helper: division algorithm on an elliptic curve
# ======================================================================

def divide_ec(P, ell, E):
    """
    Given P in E(F_q) with ell | ord(P), find Q in E(F_q) with [ell]Q = P.
    Uses exhaustive search (for small ell and small curves).
    Returns None if no such Q exists.
    """
    for Q in E.points():
        if ell * Q == P:
            return Q
    return None

# ======================================================================
# 1. Elliptic curve: 2-Sylow computation
# ======================================================================

# Use a curve over F_p where we know the group structure explicitly.
p_s = 101
F_s = GF(p_s)
lam1, lam2, lam3 = F_s(0), F_s(3), F_s(7)
E_s = EllipticCurve(F_s, [0, -(lam1+lam2+lam3), 0,
                           lam1*lam2+lam1*lam3+lam2*lam3, -lam1*lam2*lam3])
T1, T2 = E_s(lam1, 0), E_s(lam2, 0)   # basis for E[2](F_p)
ell_s = 2

def tate2_s(Ki, P):
    return int(F_s(P[0] - Ki[0])^((p_s-1)//2))

def profile2_s(P):
    return (tate2_s(T1, P), tate2_s(T2, P))

trivial_s = (1, 1)

print(f"=== Section 5.5 (EC): 2-Sylow torsion computation ===")
print(f"E: y^2 = x(x-3)(x-7) over F_{p_s},  #E = {E_s.order()}")
print()

# Apply the algorithm: start from B = {T1, T2}
B = [T1, T2]
Q_list = list(B)

print("Step 0: start from A[2] basis B = {T1, T2}")
print(f"  Q1 = T1, profile = {profile2_s(T1)}")
print(f"  Q2 = T2, profile = {profile2_s(T2)}")
print()

# Divide each Q_i by ell until the profile is non-trivial
for step in range(1, 8):
    changed = False
    for idx in range(len(Q_list)):
        Qi = Q_list[idx]
        prof = profile2_s(Qi)
        if prof == trivial_s:
            # Try to divide: find R with [2]R = Qi
            R = divide_ec(Qi, ell_s, E_s)
            if R is not None and not R.is_zero():
                Q_list[idx] = R
                changed = True
    if not changed:
        break
    print(f"Step {step}: after dividing:")
    for idx, Qi in enumerate(Q_list):
        print(f"  Q{idx+1} = {Qi},  profile = {profile2_s(Qi)},  order = {Qi.order()}")
    print()

print("Final Sylow-2 basis:")
for idx, Qi in enumerate(Q_list):
    prof = profile2_s(Qi)
    print(f"  Q{idx+1}: order={Qi.order()}, profile={prof}, non-trivial={prof != trivial_s}")

print()

# ======================================================================
# 2. Jacobian: ell = 67 Sylow torsion
#
# For the Kummer surface with p = 2^246*3*67-1:
#   J(F_p) ~ (Z_{67})^? x (other)
# The algorithm uses tate_pairing / profile to check divisibility.
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

print(f"=== Section 5.5 (Jacobian): Sylow-{ell_j} torsion computation ===")

# Step 1: Find basis B = {P1, P2} of A[ell](F_q)
# (torsion_basis does this automatically)
basis_67 = K2.torsion_basis(ell_j)
print(f"Found basis for J[{ell_j}]: {len(basis_67)} points")

# Step 2: Apply division algorithm
# Division: given P with ell | ord(P), find Q with [ell]Q = P.
# On the Kummer surface, this requires knowing the embedding field and
# using the group law.  Here we demonstrate the profile-based stopping criterion.
h_j = (p_67 + 1) // ell_j

def random_order_ell(K, h, ell):
    while True:
        P = K.random_point()
        if P.on_jacobian():
            Q = h * P
            if not Q.is_zero():
                return Q

P1 = random_order_ell(K2, h_j, ell_j)
print(f"\nSample point P1 of order {ell_j}:")
print(f"  profile t_{{[{ell_j}]}}(P1) = {P1.profile(ell_j)}")
print(f"  non-trivial = {not P1.profile(ell_j).is_trivial()}")
print()
print("A non-trivial profile confirms P1 is NOT in [ell]J(F_q),")
print("so the algorithm terminates: P1 is already a Sylow generator.")

# The full algorithm would iteratively:
#   1. Check if t_{[ell]}(Q_i) is trivial (=> can divide)
#   2. Find Q such that [ell]Q = Q_i  (division step)
#   3. Repeat until all profiles are non-trivial and independent
print()
print("Algorithm summary (Section 5.5):")
print("  1. Start with a basis B of A[ell](F_q).")
print(f"  2. For each P_i in B: while profile t_{{[{ell_j}]}}(P_i) is trivial,")
print(f"     replace P_i by a root of [{ell_j}]Q = P_i.")
print("  3. Adjust for linear dependencies in the profiles.")
print(f"  4. Output basis Q_1,...,Q_r of S_{{{ell_j},q}}(J).")
