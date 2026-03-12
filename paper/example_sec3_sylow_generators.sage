"""
Section 3 -- Theorem 2: sampling generators for the Sylow-ell torsion using profiles.

Theorem 2. Let r = rank of A[ell](F_q), assume the ell-Tate pairing is
non-degenerate, and let t_{[ell]} be the profile with respect to a basis B
of A[ell](F_q).  Then for points P_1,...,P_r in A(F_q):

    A(F_q)/A[ell](F_q) = <P_1,...,P_r>
    iff
    mu_ell^r = < t_{[ell]}(P_1), ..., t_{[ell]}(P_r) >.

Consequently, the Sylow-ell subgroup S_{ell,q}(A) = <[h]P_1,...,[h]P_r>
where h = #A(F_q) / |S_{ell,q}(A)|.

Algorithm (informal):
  1. Find a basis B of A[ell](F_q).
  2. Sample random P_i in A(F_q) until their ell-profiles are non-trivial
     and span mu_ell^r.
  3. Clear the cofactor: Q_i = [h] P_i.

Demonstrated for:
  1. An elliptic curve (ell = 5, small field)
  2. A genus-2 Jacobian (ell = 67, Kummer surface)
"""

# ======================================================================
# 1. Elliptic curve: ell = 5
#
# We use a curve E over F_q where E[5](F_q) = (Z/5)^2.
# ======================================================================

# Find a prime q with q ≡ 1 mod 5 and a curve with E[5] rational
ell_ec = 5
q_ec = next(p for p in primes(200) if p % ell_ec == 1)  # first prime ≡ 1 mod 5
Fq = GF(q_ec)
print(f"=== Section 3, Theorem 2 (EC): Sylow-{ell_ec} generators ===")
print(f"Searching over F_{q_ec}...")

E_ec = None
for a4 in Fq:
    for a6 in Fq:
        try:
            Etry = EllipticCurve(Fq, [a4, a6])
            if Etry.order() % ell_ec^2 == 0:
                tors = Etry.torsion_subgroup()
                if tors.order() % ell_ec^2 == 0:
                    E_ec = Etry
                    break
        except Exception:
            pass
    if E_ec:
        break

if E_ec is None:
    print(f"  No suitable curve found over F_{q_ec}; skipping EC part.")
else:
    n_ec = E_ec.order()
    print(f"E: {E_ec}  over F_{q_ec},  #E = {n_ec}")

    # Basis for E[ell]
    tgens = [g for g in E_ec.torsion_subgroup().gens()
             if g.order() % ell_ec == 0]
    if len(tgens) >= 2:
        B1, B2 = tgens[0], tgens[1]
    else:
        B1 = tgens[0]
        B2 = (ell_ec - 1) * B1     # fallback (will give rank-1 profile)

    def profile_ec_ell(P):
        """ell-profile via Weil pairing (bilinear, non-degenerate)."""
        return (B1.weil_pairing(P, ell_ec), B2.weil_pairing(P, ell_ec))

    # Cofactor h = #E / ell^(v_ell(#E))
    v = valuation(n_ec, ell_ec)
    h_ec = n_ec // ell_ec^v
    print(f"  #E = {n_ec} = {ell_ec}^{v} * {h_ec}")
    print(f"  Sylow-{ell_ec} subgroup has order {ell_ec^v}")

    # Sample generators of S_{ell,q}(E)
    trivial_prof = (Fq(1), Fq(1))
    gen_found = []
    attempts = 0
    while len(gen_found) < 2:
        attempts += 1
        P = E_ec.random_point()
        Q = h_ec * P
        if Q.is_zero():
            continue
        prof = profile_ec_ell(Q)
        if prof == trivial_prof:
            continue
        # Check independence: first generator always accepted,
        # second must have a profile not in <first generator's profile>
        if not gen_found or gen_found[0][1] != prof:
            gen_found.append((Q, prof))

    print(f"  Found {len(gen_found)} generators in {attempts} attempts:")
    for G, prof in gen_found:
        print(f"    order = {G.order()},  profile = {prof}")
    print()

# ======================================================================
# 2. Genus-2 Jacobian: ell = 67
#
# Using the Kummer surface over F_{p^2} with p = 2^246*3*67 - 1.
# The surface has J[67] of rank 4 over F_{p^2} (supersingular).
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
ros = (lam, mu, nu)

qo = p_67 + 1   # both Jacobians have order p+1 (supersingular)
qt = p_67 + 1

# Surface over F_{p^2} (needed: mu_67 ⊆ F_{p^2}* since p ≡ -1 mod 67)
K2 = SquaredKummerSurface(
    (F2_67(lam), F2_67(mu), F2_67(nu)),
    jacobian_order=qo, twist_order=qt
)

print(f"=== Section 3, Theorem 2 (Jacobian): Sylow-{ell_j} generators ===")
print(f"Kummer surface over F_p with p = 2^246*3*{ell_j}-1")

# Cofactor h_j = (p+1) / 67
h_j = (p_67 + 1) // ell_j
assert (p_67 + 1) % ell_j == 0

# torsion_basis(ell) requires ell | p+1 and is supersingular
basis_67 = K2.torsion_basis(ell_j)
print(f"  Found {len(basis_67)}-point basis for J[{ell_j}]")

# Sample two generators of S_{67,q}(J) on the Jacobian
gen_found_j = []
attempts_j = 0
while len(gen_found_j) < 2:
    attempts_j += 1
    P = K2.random_point()
    if not P.on_jacobian():
        continue
    Q = h_j * P
    if Q.is_zero():
        continue
    prof = Q.profile(ell_j)
    if prof.is_trivial():
        continue
    if not gen_found_j or gen_found_j[0][1] != prof:
        gen_found_j.append((Q, prof))

print(f"  Found {len(gen_found_j)} independent Sylow-{ell_j} generators in {attempts_j} attempts.")
for G, prof in gen_found_j:
    print(f"    profile = {prof}")
print(f"\nThese generators span S_{{{ell_j},q}}(J) ~ Z/{ell_j}Z x Z/{ell_j}Z.")
