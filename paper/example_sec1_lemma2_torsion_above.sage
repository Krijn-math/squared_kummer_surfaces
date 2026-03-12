"""
Section 1 -- Lemma 2 (Corte-Real Santos - Meyer - Reijnders, CEMR24):
Determining which 2-torsion point P lies above via 2-Tate pairings.

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

# ======================================================================
# 1. Elliptic curve
# ======================================================================

p = 223   # = 2^5 * 7 - 1
f = 5
h = 7

Fp  = GF(p)
Fp2.<ii> = GF(p^2)

# Re-use: find supersingular Montgomery A
A_ec = None
for A_int in range(1, p):
    A_try = Fp(A_int)
    if A_try == 2 or A_try == Fp(-2):
        continue
    E_try = EllipticCurve(Fp, [0, A_try, 0, 1, 0])
    if E_try.order() == p + 1:
        A_ec = A_try
        break

EA = EllipticCurve(Fp2, [0, Fp2(A_ec), 0, 1, 0])
print(f"Supersingular Montgomery E_A over F_p^2, A = {A_ec}, p = {p} = 2^{f}*{h}-1")

# Compute alpha: root of x^2 + Ax + 1 = 0 in F_{p^2}
R2 = Fp2['x']; xv = R2.gen()
roots_alpha = (xv^2 + A_ec*xv + 1).roots()
alpha = Fp2(roots_alpha[0][0])
print(f"alpha = {alpha}  (root of x^2 + {A_ec}x + 1 = 0 in F_p^2)")

T0 = EA(Fp2(0), Fp2(0))
T1 = EA(alpha, 0)
T2 = EA(1/alpha, 0)
print(f"T0 = (0,0),  T1 = (alpha, 0),  T2 = (1/alpha, 0)")

def tate2_fp2(Ti, P):
    """
    2-Tate pairing t_2(T_i, P) over F_{p^2}.
    For the reduced Tate pairing of degree 2 over F_{p^2}:
      t_2(T_i, P) = (x_P - x_{T_i})^{(p^2-1)/2}  in {+1,-1}.
    """
    diff = P[0] - Ti[0]
    return int(Fp2(diff)^((p^2 - 1) // 2))

print("\nVerifying Lemma 2: [2^{f-1}]P = T_i  iff  t_2(T_i,P)=1, t_2(T_j,P)=-1 for j!=i")
print()

# Sample points of order 2^f in E(F_{p^2}) and check which T_i they are above
# First, get a basis for E[2^f]: clear the h-cofactor from a random point
def random_2f_point(E, f, h):
    """Return a random point of order exactly 2^f."""
    while True:
        Q = E.random_point()
        Q = h * Q      # clear cofactor; now order divides 2^f
        if Q.order() == 2^f:
            return Q

for trial in range(6):
    P = random_2f_point(EA, f, h)
    half_P = (2^(f-1)) * P        # = T_i for some i

    # Determine which T_i
    if half_P == T0:
        expected_above = 0
    elif half_P == T1:
        expected_above = 1
    elif half_P == T2:
        expected_above = 2
    else:
        expected_above = -1        # T0, T1, T2 are defined over F_{p^2}

    t_vals = (tate2_fp2(T0, P), tate2_fp2(T1, P), tate2_fp2(T2, P))
    pairing_says = t_vals.index(1)   # the unique T_i with t_2(T_i, P) = 1

    match = (pairing_says == expected_above)
    print(f"  [2^{f-1}]P = T{expected_above},  t_2 values = {t_vals},  "
          f"pairing says T{pairing_says}  {'OK' if match else 'FAIL'}")

print()

# ======================================================================
# 2. Jacobian analogue (Gaudry-Schost)
#
# J[2] has 4 basis elements.  The 2-Tate profile with respect to the
# 4-element basis encodes which 2-torsion point [2^{k-1}]P equals:
# a unique basis element has pairing value True (square), all others False.
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
basis_gs = K_gs.two_torsion_basis()  # T_12, T_13, T_14, T_15

print("=== Jacobian (Gaudry-Schost) analogue of Lemma 2 ===")
print("2-profile of P  =>  which 2-torsion point [2^{k-1}]P equals.")
print()

trivial = [True, True, True, True]
for _ in range(5):
    P  = K_gs.random_point()
    P2 = 2 * P
    P4 = 2 * P2  # = 4P  -- moving toward [2^{k-1}]P
    prof_P  = P.two_profile(basis_gs)
    prof_2P = P2.two_profile(basis_gs)

    # Profile of P tells us 4 independent pairing values; comparing P and 2P:
    print(f"  P  profile: {prof_P}")
    print(f"  2P profile: {prof_2P}  (trivial={prof_2P == trivial})")
    print()
