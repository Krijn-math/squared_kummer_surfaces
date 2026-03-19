"""
In this example, we showcase the entangled basis generation technique by [ZSDPB18],
although the interpretation in terms of profiles was clarified later by Damien Robert.

At first, this technique looks magical, but once we look at the 2-profiles, it suddenly
becomes clear! There are a few things important though: First, we work with a supersingular
curve, so that we have symmtric Sylow-2 torsion! Second, we work with Montgomery curves,
although so that T_0 = (0,0) is a 2-torsion point.

Both statements may be generalized, e.g., we need a curve with symmetric Sylow-2 torsion
and any curve model where (0,0) is a 2-torsion point, and generalized even more as we 
show later.

Section 1 -- Lemma 1 (Zanon-Simplicio-Pereira-Doliskani-Barreto, ZSPDB18):
Entangled basis for the 2^f-torsion of a maximal supersingular Montgomery curve.

Lemma 1. Let E_A : y^2 = x^3 + Ax^2 + x be a maximal supersingular Montgomery
curve over F_{p^2} with p = 2^f * h - 1 (f in N, h odd cofactor, A != 0).
Let t in F_{p^2} be a non-square such that x_P = -A/(1+t^2) is also non-square.
Then:
  - x_P defines a point P in E(F_{p^2}),
  - x_Q = -x_P - A defines a point Q in E(F_{p^2}),
  - ([h]P, [h]Q) is a basis for E[2^f].

The key insight is that choosing x_P as a non-square ensures P has a non-trivial
2-profile, and the specific choice x_Q = -x_P - A gives Q a *different* non-trivial
profile, so the two points span different cosets of E/[2]E.
"""

# ======================================================================
# 1. Elliptic curve: supersingular Montgomery curve over F_{p^2}
# ======================================================================
# Use p = 223 = 2^5 * 7 - 1, so f = 5, h = 7.
# Find a supersingular Montgomery curve over F_p (trace = 0 => #E = p+1 = 224).
print()
print("=== Section 1: The magic of entangled basis generation")
print()

p = 223   # = 2^5 * 7 - 1
f = 5
h = 7
assert p == 2^f * h - 1 and is_prime(p)

Fp  = GF(p)
R_p = Fp['x']; x = R_p.gen()

# Search for supersingular Montgomery A over F_p (need #E_A(F_p) = p+1)
A = 6
print(f"Supersingular Montgomery curve: A = {A},  p = {p} = 2^{f}*{h}-1")

# Work over F_{p^2} where E[2^f] is fully rational
Fp2.<ii> = GF(p^2)
EA = EllipticCurve(Fp2, [0, Fp2(A), 0, 1, 0])

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

# first, we show that any point with non-square x-coordinate has a non-trivial profile
# which is obvious when you consider that the squareness of this x-coordinate is 
# exactly the 2-Tate pairing with T0 = (0,0)

# we take a global non-square z ∈ Fp2 and then sample non-squares below by z*a^2 for any a ∈ Fp2
while True:
    Z = Fp2.random_element()
    if not Z.is_square():
        break

for i in range(100):
    while True:
        a = Fp2.random_element()
        xp = Z*a**2
        if (xp**3 + A*xp**2 + xp).is_square():
            P = EA.lift_x(xp)
            break
    
    assert two_profile(P) != (1,1,1)

    # in fact, we get that the order of P is divisible by 2^f, and doesn't lie above T0
    Q = ((p+1) // 2)*P
    assert Q == T1 or  Q == T2


# now we demonstrate the magic:
# if we sample xp in the right way, we can define an xq that defines a point too, with a 
# different non-trivial profile!

# Find a non-square t in F_{p^2} such that x_P = -A/(1+t^2) is also non-square in F_{p^2}
# and defines a point on E, then define P and Q

while True:
    t = Fp2.random_element()
    t = Z*t**2
    xP = -A / (1 + t**2)            # in practice, you can precompute a list of t to ensure xP is non-square
    if not xP.is_square() and (xP**3 + A*xP**2 + xP).is_square():
        break

xQ = -xP - A
P = EA.lift_x(xP)
Q = EA.lift_x(xQ)

print(f"\nChosen non-square t   =   {t}")
print(f"x_P =   -A/(1+t^2)    =   {xP}")
print(f"x_Q =   -x_P - A      =   {xQ}")
print()

print(f"profile of P: {two_profile(P)}")
print(f"profile of Q: {two_profile(Q)}")

# Now we show that these indeed generate the 2^f-torsion
# Clear the h-cofactor to land on the 2^f-torsion
P = h * P
Q = h * Q

assert P.order() == 2**f and Q.order() == 2**f, "Points should have order 2^f"
print(f"\n[h]P order = {factor(P.order())}")
print(f"[h]Q order = {factor(Q.order())}")
print()


# To verify ([h]P, [h]Q) is a basis for E[2^f], there are ample things we can do
# the fact that their profiles are different essentially shows it already
# but lets do some more

# Check: [2^(f-1)] P and [2^(f-1)] Q are distinct two-torsion points
assert (2**(f-1)) * P != (2**(f-1)) * Q

# Check: Weil pairing has order 2^5
zeta = P.weil_pairing(Q, 2**f)
print(f"The Weil pairing e_2^f(P, Q) has order {factor(zeta.multiplicative_order())}")


# ======================================================================
#       FOR HIGHER DIMENSIONS
# Surprisingly, this trick only seems to work for ell = 2 and dim = 1!
# I don't know of any work that generalizes this to ell > 2 or higher
# dimensions. Future work may focus on the specific curve model that
# enables this trick on elliptic curve. So, for ell = 3, one may consider
# Hessians, and similarly for higher dimensions, so that one can easily
# describe at least one Tate pairing of degree ell almost trivially.
# ======================================================================


"""
Given what we know, we may now generalize the above technique to any curve
model where the 2-Tate pairings are "easy", although we need q to be even
so that -1 is a square.

Section 5.6 -- Generalized entangled basis generation (Lemma 5).

Lemma 5. Let E : y^2 = (x - lam1)(x - lam2)(x - lam3) with lam_i in F_q,
and let h = #E(F_q) / 2^{f+g}  (cofactor).  For a non-square u in F_q, set

    x_P = (lam2 + lam3 - 2*lam1) / (1 + u^2).

If x_P is a non-square and defines a point P in E(F_q), then
    x_Q = -x_P + lam2 + lam3
defines a point Q in E(F_q) such that [h]P and [h]Q generate the
Sylow-2 subgroup  S_{2,F_q}(E) ≅ Z_{2^f} x Z_{2^g}  with f >= g.

Key identities:
    x_Q - lam2 = -(x_P - lam3),      x_Q - lam3 = -(x_P - lam2).
These permute the quadratic characters, ensuring t_{[2]}(P) != t_{[2]}(Q).

This generalizes Lemma 1: the Montgomery case uses lam1 = 0 and
    x_P = -A/(1+t^2) = -(lam2+lam3)/(1+t^2) = (lam2+lam3-2*0)/(1+t^2).
"""

print()    
print()    
print()
print("=== Section 5: Generalizing entangled basis generation")
print()

# ======================================================================
# 1. Elliptic curve: Lemma 5 construction
# ======================================================================

# Use a concrete curve with three rational 2-torsion points.
# Pick lam1=0, lam2=3, lam3=7 over F_101 (from Theorem 1 example).
p = 4999
F = GF(p)
F2 = GF(p**2)

# we take any random curve, and check what the Sylow-2 torsion is
lam1 = F.random_element()
lam2 = F.random_element()
lam3 = F.random_element()

E = EllipticCurve(F2, [0, -(lam1+lam2+lam3), 0,
                    lam1*lam2+lam1*lam3+lam2*lam3, -lam1*lam2*lam3])

# here, h is the cofactor which kills everything but the 2^k torsion
h = E.order() // 2**(valuation(E.order(), 2))
P, Q = E.abelian_group().gens()
P = h*P
Q = h*Q
first_sylow = P.order()
second_sylow = Q.order()
print(f"E: y^2 = (x-{lam1})*(x-{lam2})*(x-{lam3}) over F_{p}^2")
print(f"E(Fq) has Sylow-2 torsion {factor(first_sylow)} x {factor(second_sylow)}")
print()    

# of course we know the 2-torsion and their pairings
T1 = E(lam1, 0)
T2 = E(lam2, 0)
T3 = E(lam3, 0)

def two_profile(P):
    # Tate profile of degree 2, e.g. all tate pairings of P with respect to the 2-torsion on E
    return ( T1.tate_pairing(P, 2, 2), T2.tate_pairing(P, 2, 2), T3.tate_pairing(P, 2, 2) )   

# first take any non-square
while True:
    Z = F2.random_element()
    if not Z.is_square():
        break

# Then, find u in F_p2 such that x_P - lam1 = (lam2 + lam3 - 2*lam1) / (1 + u^2)
# is non-square
delta = (lam2 + lam3 - 2*lam1)

while True:
    u = F2.random_element()
    t = delta * u**2 / (1 + u**2)
    t = delta / (1 + u**2)
    xP = t + lam1
    # we want xP - lam1 non-square, that is,  delta / (1 + u**2) is non-square
    # and furthermore that xP defines a point on the curve
    if not (delta/(1+u**2)).is_square() and ((xP - lam1)*(xP - lam2)*(xP - lam3)).is_square():
        break

# now, "magically", xQ also defines a point and P, Q will generate the Sylow-2 torsion
xQ = -xP + lam2 + lam3
P = E.lift_x(xP)
Q = E.lift_x(xQ)


print(f"\nChosen u   =   {u}")
print(f"x_P = (lam2 + lam3 - 2*lam1) / (1 + u^2) + lam1 = {xP}")
print(f"x_Q =                         -xP + lam2 + lam3 = {xQ}")
print()

prof_P = two_profile(P)
prof_Q = two_profile(Q)
assert prof_P != (1,1,1)
assert prof_Q != (1,1,1)
assert prof_P != prof_Q

print(f"2-profile of P: {prof_P}")
print(f"2-profile of Q: {prof_Q}")
print()

# clear the cofactor to get a basis of the Sylow-2 torsion
P = h * P
Q = h * Q

# however, very likely, these will not be presicely of the orders
# so we can clean up the representation a bit by a trick

# first, we swap to ensure P has order first_sylow
if P.order() < Q.order():
    R = P
    P = Q
    Q = R
    
assert P.order() == first_sylow 

# then, adjust Q to have order second_sylow, but still non-trivial profile!
R = Q
while True:
    R += P
    if R.order() == second_sylow and not two_profile(R) == (1,1,1):
        Q = R
        break

print(f"\n[h]P order = {factor(P.order())}")
print(f"[h]Q order = {factor(Q.order())}")
print()

# now we have a nice basis for our sylow-2 torsion using the generalied entangled basis generation!

print("Verified: the 'magic' of entangled basis generation is just 2-profiles!")
print("Furthermore, if we know enough about the 2-torsion, we can generalize")
print("this technique beyond supersingular curves and Montgomery curves!\n")