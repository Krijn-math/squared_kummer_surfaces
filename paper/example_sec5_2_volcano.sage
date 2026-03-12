"""
Section 5.2 + Example 3 -- Pairing the volcano.

For a prime p ≡ 7 mod 8, supersingular Montgomery curves E_A : y^2 = x^3+Ax^2+x
over F_p form a 2-volcano:
  - Surface: curves with End_p(E) ≅ Z[(1+sqrt(-p))/2]  (E[2](F_p) = Z/2 x Z/2)
  - Floor:   curves with End_p(E) ≅ Z[sqrt(-p)]          (E[2](F_p) = Z/2)

The self-Tate pairing t_2(T_i, T_i) distinguishes the isogeny type of T_i:
  t_2(T_i, T_i) = -1   =>  T_i generates a *vertical* (descending) isogeny
  t_2(T_i, T_i) =  1   =>  T_i generates a *horizontal* isogeny

Self-pairing formula (via bilinearity):
  t_2(T_i, T_i) = t_2(T_i, T_i + T_j) / t_2(T_i, T_j)
                = Legendre(x_{T_i+T_j} - x_{T_i}) / Legendre(x_{T_j} - x_{T_i})

where T_i + T_j is the third 2-torsion point.

Example 3 from the paper: p = 23, A = 10.
"""

# ======================================================================
# Example 3 from the paper: p = 23, A = 10
# ======================================================================

p = 23
A = 10
F = GF(p)

E10 = EllipticCurve(F, [0, A, 0, 1, 0])   # E_A : y^2 = x^3 + 10x^2 + x
print(f"=== Section 5.2, Example 3: pairing the volcano ===")
print(f"p = {p}  (≡ {p % 8} mod 8),  A = {A}")
print(f"Curve E_{A}: y^2 = x^3 + {A}x^2 + x over F_{p}")
print(f"#E_{A}(F_{p}) = {E10.order()}")
print()

# 2-torsion points: y = 0, so x^3 + Ax^2 + x = x(x^2+Ax+1) = 0
# => T0 = (0,0) and roots of x^2 + Ax + 1 = 0
T0 = E10(0, 0)
R = F['x']; xv = R.gen()
roots_other = (xv^2 + A*xv + 1).roots()
assert len(roots_other) == 2, f"Both roots should be in F_{p}"
T17, T19 = E10(roots_other[0][0], 0), E10(roots_other[1][0], 0)
# Sort so T17 has x=17, T19 has x=19
if T17[0] > T19[0]:
    T17, T19 = T19, T17
print(f"2-torsion points: T0 = {T0},  T{T17[0]} = {T17},  T{T19[0]} = {T19}")
print()

# 2-torsion group: T0 + T17 + T19 = O => T0+T17 = T19, T0+T19 = T17, T17+T19 = T0
assert T0 + T17 + T19 == E10(0)

def legendre(x, p):
    return int(GF(p)(x)^((p-1)//2))

def tate2_pair(Ti, Tj, p):
    """t_2(T_i, T_j) = Legendre(x_{T_j} - x_{T_i}, p) for i != j."""
    return legendre(Tj[0] - Ti[0], p)

def self_tate2(Ti, Tj, p):
    """
    Self-pairing t_2(T_i, T_i) via bilinearity:
      = t_2(T_i, T_i + T_j) / t_2(T_i, T_j)
    where T_i + T_j is the third 2-torsion point.
    """
    Ti_plus_Tj = Ti + Tj   # = the third 2-torsion point
    num = tate2_pair(Ti, Ti_plus_Tj, p)
    den = tate2_pair(Ti, Tj, p)
    return num * den   # in {+1,-1}: division is the same as multiplication

tors = [T0, T17, T19]
names = {T0: "T0", T17: f"T{T17[0]}", T19: f"T{T19[0]}"}

print("Self-pairings t_2(T_i, T_i):")
for Ti in tors:
    # Use any other T_j != T_i for the formula
    Tj = next(t for t in tors if t != Ti)
    sp = self_tate2(Ti, Tj, p)
    isogeny_type = "VERTICAL (descending)" if sp == -1 else "horizontal"
    print(f"  t_2({names[Ti]}, {names[Ti]}) = {sp}  =>  {names[Ti]} generates a {isogeny_type} isogeny")

print()
print("Verification against paper:")
print("  T0 and T17 should give self-pairing +1 (horizontal),")
print("  T19 should give self-pairing -1 (vertical).")
assert self_tate2(T0, T17, p) == 1
assert self_tate2(T17, T0, p) == 1
assert self_tate2(T19, T0, p) == -1
print("  All checks passed.")
print()

# ======================================================================
# Compute the three 2-isogenies from E10 and verify the codomain curves
# ======================================================================

print("The three 2-isogenies from E_10:")
isogens = E10.isogenies_prime_degree(2)
for phi in isogens:
    cod = phi.codomain()
    ker_gen = phi.kernel_polynomial().roots(F)
    # The kernel generator is the 2-torsion point with that x-coordinate
    ker_x = ker_gen[0][0] if ker_gen else None
    # Find which T_i this corresponds to
    ker_name = next((names[t] for t in tors if t[0] == ker_x), f"x={ker_x}")
    # Get A-coefficient of codomain
    cod_A = cod.a2()
    print(f"  phi_{ker_name}: E_{A} -> E_{cod_A}  (kernel generator = {ker_name})")
    if ker_x == T19[0]:
        print(f"    -> vertical isogeny (T19 has self-pairing = -1)")
    else:
        print(f"    -> horizontal isogeny")

print()

# ======================================================================
# General lesson: for p ≡ 7 mod 8
# T_i generates a vertical isogeny  iff  t_2(T_i, T_i) = -1
# This follows because the vertical isogeny's kernel point is above
# the *non-square* directions, while horizontal isogenies preserve the
# surface structure.
# ======================================================================

# Also show using a general p ≡ 7 mod 8 example
for p2 in primes(100):
    if p2 % 8 != 7:
        continue
    F2 = GF(p2)
    for A2_int in range(1, p2):
        A2 = F2(A2_int)
        if A2 in (F2(2), F2(-2)):
            continue
        E2 = EllipticCurve(F2, [0, A2, 0, 1, 0])
        if E2.order() != p2 + 1:
            continue
        R2 = F2['x']; xv2 = R2.gen()
        rt = (xv2^2 + A2*xv2 + 1).roots()
        if len(rt) < 2:
            continue
        T0_2 = E2(0, 0)
        Ta_2 = E2(rt[0][0], 0)
        Tb_2 = E2(rt[1][0], 0)
        sp_a = self_tate2(Ta_2, Tb_2, p2)
        sp_b = self_tate2(Tb_2, Ta_2, p2)
        # Exactly one should be -1 (the vertical isogeny)
        assert sorted([sp_a, sp_b]) == [-1, 1], f"Expected exactly one vertical at p={p2}"
        break
    break
print(f"Spot check for p={p2}, A={A2}: exactly one 2-torsion point generates a vertical isogeny. OK")
