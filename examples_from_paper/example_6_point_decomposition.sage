"""
Section 5.1 -- Decomposing R as aP + bQ via Tate pairings.

Let E/F_{p^2} be a maximal supersingular elliptic curve with 2^f | p+1,
and let (P, Q) be a basis for E[2^f].  Given a third point R in E[2^f],
we want to find a, b such that R = aP + bQ.

Step 1: Compute profiles
    t_{[2^f]}(P) = (t_{2^f}(P,P), t_{2^f}(Q,P))
    t_{[2^f]}(Q) = (t_{2^f}(P,Q), t_{2^f}(Q,Q))
    t_{[2^f]}(R) = (t_{2^f}(P,R), t_{2^f}(Q,R))

Step 2: Use bilinearity t_{2^f}(K, aP+bQ) = t_{2^f}(K,P)^a * t_{2^f}(K,Q)^b
    to solve for a, b as discrete logarithms:
    a = dlog_{zeta_PP}(t_{2^f}(P,R))   where zeta_PP = t_{2^f}(P,P)
    b = dlog_{zeta_QQ}(t_{2^f}(Q,R))   (using independent second component)

Demonstrated for:
  1. Elliptic curve (using the Weil pairing as proxy for the Tate pairing)
  2. Genus-2 Jacobian (ell = 67 Kummer surface)
"""

# ======================================================================
# 1. Elliptic curve: decompose R = aP + bQ in E[ell]
# ======================================================================

p_ec  = 223   # = 2^5 * 7 - 1; supersingular Montgomery
f_ec  = 5
h_ec  = 7
ell_ec = 2^f_ec   # = 32

Fp  = GF(p_ec)
Fp2.<ii> = GF(p_ec^2)

# Find supersingular A
A_ec = None
for A_int in range(1, p_ec):
    A_try = Fp(A_int)
    if A_try in (Fp(2), Fp(-2)):
        continue
    if EllipticCurve(Fp, [0, A_try, 0, 1, 0]).order() == p_ec + 1:
        A_ec = A_try
        break

EA = EllipticCurve(Fp2, [0, Fp2(A_ec), 0, 1, 0])

# Basis (P, Q) for E[2^f] via Lemma 1 construction
def legendre_fp(x):
    return int(Fp(x)^((p_ec-1)//2))

for t_int in range(1, p_ec):
    t = Fp(t_int)
    if legendre_fp(t) == -1:
        xP = -A_ec / (1 + t^2)
        if xP != 0 and legendre_fp(xP) == -1:
            break

xQ_coord = -xP - A_ec
P_ec = EA.lift_x(Fp2(xP))
Q_ec = EA.lift_x(Fp2(xQ_coord))
P_ec = h_ec * P_ec
Q_ec = h_ec * Q_ec

assert P_ec.order() == ell_ec and Q_ec.order() == ell_ec, "Basis points should have order 2^f"

print(f"=== Section 5.1 (EC): decompose R = aP + bQ in E[2^{f_ec}] ===")
print(f"Supersingular Montgomery, p = {p_ec}, f = {f_ec}")

# Weil pairing as proxy for Tate pairing (bilinear, non-degenerate over F_{p^2})
zeta_PP = P_ec.weil_pairing(P_ec, ell_ec)   # Note: weil_pairing is alternating so this = 1
# For weil_pairing, e(P,P) = 1 (alternating).  We use e(P,Q) as the generator.
zeta_PQ = P_ec.weil_pairing(Q_ec, ell_ec)
print(f"  zeta = e(P, Q) = {zeta_PQ}  (should be a primitive {ell_ec}-th root of unity)")
print(f"  zeta^{ell_ec} = {zeta_PQ^ell_ec}  (should be 1)")

# Choose random a, b and form R = aP + bQ
a_true = randint(1, ell_ec - 1)
b_true = randint(1, ell_ec - 1)
R_ec = a_true * P_ec + b_true * Q_ec

# Decompose R: use t(P, R) = e(P,Q)^b and t(Q, R) = e(P,Q)^{-a}
# (since weil is alternating: e(P, aP+bQ) = e(P,Q)^b)
t_PR = P_ec.weil_pairing(R_ec, ell_ec)   # = zeta^b
t_QR = Q_ec.weil_pairing(R_ec, ell_ec)   # = zeta^{-a}

b_rec = discrete_log(t_PR, zeta_PQ, ell_ec)
a_rec = discrete_log(t_QR^(-1), zeta_PQ, ell_ec)

print(f"\n  True  (a, b) = ({a_true}, {b_true})")
print(f"  Recovered (a, b) = ({a_rec}, {b_rec})")
assert a_true % ell_ec == a_rec % ell_ec and b_true % ell_ec == b_rec % ell_ec
print(f"  Verification: {a_rec}*P + {b_rec}*Q == R: {a_rec*P_ec + b_rec*Q_ec == R_ec}")
print()

# ======================================================================
# 2. Jacobian: decompose R = aP + bQ in J[ell] with ell = 67
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

print(f"=== Section 5.1 (Jacobian): decompose R = aP + bQ in J[{ell_j}] ===")

# Get two independent order-67 points on the Jacobian
h_j = (p_67 + 1) // ell_j

def random_order_ell(K, h, ell):
    while True:
        P = K.random_point()
        if not P.on_jacobian():
            continue
        Q = h * P
        if not Q.is_zero():
            return Q

P_j = random_order_ell(K2, h_j, ell_j)
Q_j = random_order_ell(K2, h_j, ell_j)
# Ensure independence: Q_j not a multiple of P_j
while True:
    zeta_test = P_j.tate_pairing(Q_j, ell_j)
    if zeta_test != 1:
        break
    Q_j = random_order_ell(K2, h_j, ell_j)

zeta_j = P_j.tate_pairing(Q_j, ell_j)
print(f"  zeta = t_{ell_j}(P, Q) has order {zeta_j.multiplicative_order()}")

# Choose random a, b and form R = aP + bQ
a_j = randint(1, ell_j - 1)
b_j = randint(1, ell_j - 1)
R_j = a_j * P_j + b_j * Q_j

# Recover a, b via pairings
# t(P, R) = t(P, aP+bQ) = t(P,P)^a * t(P,Q)^b = zeta^b  (assuming t(P,P)=1 on Kummer)
# t(Q, R) = t(Q, aP+bQ) = t(Q,P)^a * t(Q,Q)^b = zeta^{-a}
t_PR_j = P_j.tate_pairing(R_j, ell_j)
t_QR_j = Q_j.tate_pairing(R_j, ell_j)

b_j_rec = discrete_log(t_PR_j, zeta_j, ell_j)
# t(Q,P) = zeta^{-1}  (Weil antisymmetry inherited by Tate)
t_QP_j = Q_j.tate_pairing(P_j, ell_j)
a_j_rec = discrete_log(t_QR_j, t_QP_j, ell_j)

print(f"\n  True  (a, b) = ({a_j}, {b_j})")
print(f"  Recovered (a, b) = ({a_j_rec}, {b_j_rec})")
check = (a_j_rec * P_j + b_j_rec * Q_j == R_j)
print(f"  Verification: {a_j_rec}*P + {b_j_rec}*Q == R: {check}")
