p = 13
F.<z2> = GF(p**2)
R.<x> = F[]

λ = 8*z2 + 5
μ = 5*z2
ν = 5*z2 + 5

# while True:
#     used = [0,1]
    
#     while True:
#         λ = F.random_element()
#         if not λ in used:
#             used.append(λ)
#             break
    
#     while True:
#         μ = F.random_element()
#         if not μ in used:
#             used.append(λ)
#             break
    
#     while True:
#         ν = F.random_element()
#         if not ν in used:
#             used.append(λ)
#             break        

H = HyperellipticCurve(x*(x-1)*(x-λ)*(x-μ)*(x-ν))
J = H.jacobian()

def random_curve_point(curve):
    while True:
        try:
            return curve.lift_x(F.random_element())
        except:
            pass
        
def random_jacobian_element(jac):
    P1 = random_curve_point(jac.curve())
    P2 = random_curve_point(jac.curve())
    return jac(P1) + jac(P2)

def major_order(P):
    zero = P.scheme()(0)
    n = 0
    while True:
        n += 1
        if factorial(n)*P == zero:
            return n

def max_val(P):
    full = prod(primes(2, major_order(P)))
    k = 0
    while P != 0*P:
        k += 1
        P = full*P
    
    return k   

def order(P):
    zero = P.scheme()(0)    
    ord = 1
    Q = P
    while Q != zero:
        ord +=1
        Q = Q + P
        
    return ord

def sample_order(ell):
    for i in range(10):
        P = random_jacobian_element(J)
        ord = order(P)
        if ord % ell == 0:
            return (ord // ell)*P

points = []
orders = []
prime_facs = []
for i in range(10):
    P = random_jacobian_element(J)
    z = order(P)
    
    points.append(P)
    orders.append(z)
    prime_facs += prime_factors(z)

prime_facs = Set(prime_facs).list()
# prime_facs.pop(prime_facs.index(2))
max_vals = [ max([ valuation(ord, ell) for ord in orders ]) for ell in prime_facs ]

print(orders)
print(prime_facs)
print(max_vals)
print()

    # high = max(max_vals)
    # if high > 2:
    #     break

assert 3 == 5


"""

TODO: the dim-2 example should really be on a more interesting surface, we can take anything essentially

Section 2.2 -- Profiles assuming rational n-torsion.

When E[n](F_q) = E[n]  (full n-torsion is F_q-rational), the profile map

    t_{[n]} : E(F_q)/[n]E(F_q)  ~->  mu_n^r

becomes an *isomorphism*.  With a basis (K_1,...,K_r) of E[n](F_q), the
profile completely identifies the coset of any P in E(F_q).

For n=2 with r=2 (maximal): the profile t_{[2]}(P) = (t_2(K_1,P), t_2(K_2,P))
gives an isomorphism E(F_q)/[2]E(F_q) ~-> {+1,-1}^2 ~ Z/2 x Z/2.

The third 2-Tate pairing t_2(K_3, P) (for K_3 = K_1 + K_2) is determined
by the first two since t_2 is bilinear:
    t_2(K_1 + K_2, P) = t_2(K_1, P) * t_2(K_2, P).

Demonstrated for:
  1. An elliptic curve with full rational 2-torsion and full rational 5-torsion
  2. A genus-2 Jacobian (Gaudry-Schost Kummer surface, 2-torsion case)
"""

# ======================================================================
# 1a. Elliptic curve: degree n=2 with E[2] fully rational
# ======================================================================

p = 101
F = GF(p)
lam1, lam2, lam3 = F(0), F(3), F(7)
E = EllipticCurve(F, [0, -(lam1+lam2+lam3), 0,
                      lam1*lam2 + lam1*lam3 + lam2*lam3, -lam1*lam2*lam3])
K1, K2 = E(lam1, 0), E(lam2, 0)   # basis for E[2](F_p) ~ Z/2 x Z/2
K3 = K1 + K2                        # = (lam3, 0) = third 2-torsion point

def tate2(Ki, P, p):
    return int(GF(p)(P[0] - Ki[0])^((p-1)//2))

print("=== Section 2.2 (EC, n=2): profile isomorphism E(F_p)/[2]E(F_p) ~-> {+1,-1}^2 ===")
print(f"E: y^2 = x(x-3)(x-7) over F_{p},  E[2](F_p) = Z/2 x Z/2")
print(f"Basis K1=(0,0), K2=(3,0);  K3=K1+K2={K3}")
print()

# All four cosets, each identified by a unique profile
cosets = {}
for pt in E.points():
    if pt.is_zero():
        continue
    prof = (tate2(K1, pt, p), tate2(K2, pt, p))
    cosets.setdefault(prof, []).append(pt)

print("Profile isomorphism (2 independent pairings suffice):")
for prof, pts in sorted(cosets.items()):
    print(f"  t_{{[2]}}(P) = {prof}: {len(pts)} points")

print()
# Verify bilinearity: t_2(K3, P) = t_2(K1,P)*t_2(K2,P)
for pt in E.points():
    if pt.is_zero():
        continue
    lhs = tate2(K3, pt, p)
    rhs = tate2(K1, pt, p) * tate2(K2, pt, p)
    assert lhs == rhs, f"Bilinearity failed at {pt}"
print("Verified: t_2(K1+K2, P) = t_2(K1,P) * t_2(K2,P) for all P in E(F_p).")
print()

# ======================================================================
# 1b. Elliptic curve: degree n=5 with E[5] fully rational over some F_q
#
# We find a curve over F_q with E[5](F_q) = (Z/5Z)^2 and demonstrate
# the profile isomorphism E(F_q)/[5]E(F_q) ~-> mu_5^2.
# ======================================================================

# Need q ≡ 1 mod 5 so that mu_5 ⊆ F_q, and a curve with E[5](F_q) = (Z/5)^2.
# The smallest such prime: q = 11 (11 ≡ 1 mod 5).
# Search for E over F_11 with 5^2 | #E(F_11) (so E[5] is fully rational).

q5 = None
E5 = None
for q_try in primes(100):
    if q_try % 5 != 1:
        continue
    Fq = GF(q_try)
    for a4 in Fq:
        for a6 in Fq:
            try:
                Etry = EllipticCurve(Fq, [a4, a6])
                ord_ = Etry.order()
                if ord_ % 25 == 0 and Etry.torsion_subgroup().order() % 25 == 0:
                    q5 = q_try
                    E5 = Etry
                    break
            except Exception:
                pass
        if E5:
            break
    if E5:
        break

if E5 is not None:
    Fq = GF(q5)
    print(f"=== Section 2.2 (EC, n=5): profile isomorphism E(F_{q5})/[5]E(F_{q5}) ~-> mu_5^2 ===")
    print(f"E: {E5}  over F_{q5},  #E = {E5.order()}")
    tors5 = [pt for pt in E5.torsion_subgroup().gens()]
    print(f"5-torsion generators: {tors5}")

    # Profile with respect to a basis (K1, K2) of E[5](F_q)
    K1_5, K2_5 = tors5[0], tors5[1] if len(tors5) >= 2 else (tors5[0], 3*tors5[0])
    def profile5(P):
        return (P.weil_pairing(K1_5, 5), P.weil_pairing(K2_5, 5))

    seen5 = set()
    for pt in E5.points():
        prof = profile5(pt)
        seen5.add(prof)
    print(f"Distinct 5-profiles seen: {len(seen5)} (expected <= 25)")
    print()

# ======================================================================
# 2. Jacobian (Gaudry-Schost): n=2, rank-4 case
#
# J[2](F_p) = Z_2^4, giving 2^4=16 cosets of J(F_p)/[2]J(F_p).
# The profile map t_{[2]} : J(F_p) -> {True,False}^4 is an isomorphism.
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

print("=== Section 2.2 (Jacobian): profile isomorphism J(F_p)/[2]J(F_p) ~-> {0,1}^4 ===")
print("Gaudry-Schost: J[2](F_p) = Z_2^4,  16 cosets")
print()

trivial = [True, True, True, True]
profiles_seen = {}
n_samples = 40
for _ in range(n_samples):
    P = K_gs.random_point()
    prof = tuple(P.two_profile(basis_gs))
    profiles_seen[prof] = profiles_seen.get(prof, 0) + 1

print(f"After {n_samples} samples: {len(profiles_seen)} distinct cosets encountered.")
print("Sample of profiles:")
for prof, cnt in list(profiles_seen.items())[:6]:
    tag = " (trivial=[2]J)" if list(prof) == trivial else ""
    print(f"  {prof}: {cnt} times{tag}")

print()
print("The full isomorphism t_{[2]} : J(F_p)/[2]J(F_p) ~-> Z_2^4 would require")
print("sampling all 16 cosets, which happens quickly for random points.")
