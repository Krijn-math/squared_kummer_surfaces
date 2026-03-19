"""
SquaredKummerPoint: a point on a SquaredKummerSurface.

Supports scalar multiplication via a Montgomery-like ladder (xDBL / xDBLADD),
membership testing, normalization, and computation of the pairing profile.
"""

from operator import index as _int
from math import prod
from sage.all import factor

class TateProfile:
    """
        helper class to work with Tate profiles, as r-tuples in mu_d roots of unity
    """

    def __init__(self, coords):
        self.coords = tuple(coords)
        self.field = coords[0].base_ring()
        self.length = len(coords)
        
    def is_trivial(self):
        if self.coords == tuple([1 for i in range(self.length)]):
            return True
        
        return False
    
    def __repr__(self):
        return f"{self.coords}"
    
    def __getitem__(self, key):
        return self.coords[key]
    
    def __eq__(self, other):
        return self.coords == other.coords
        
    def __pow__(self, other):
        return TateProfile([c**other for c in self.coords])

class SquaredKummerPoint:

    def __init__(self, surface, coords, addition_matrix=None):
        """
        surface:          the SquaredKummerSurface this point lives on.
        coords:           4-tuple of field elements (X : Y : Z : T).
        addition_matrix:  optional 4x4 matrix; set when this point is a
                          2-torsion point used as a pairing input.
        """
        self.surface = surface
        self.coords = tuple(coords)
        assert self.on_kummer()
        self._addition_matrix = addition_matrix

    def _check_same_surface(self, other, name="other"):
        if not isinstance(other, SquaredKummerPoint):
            raise TypeError(f"{name} must be a SquaredKummerPoint, got {type(other)}")
        if other.surface != self.surface:
            raise ValueError(f"{name} must lie on the same SquaredKummerSurface")
        
    def same_jacobian(self, other):
        """
        Return True if self and other originate from the same Jacobian (not one
        from the Jacobian and the other from its quadratic twist).

        Detection: the coordinate X1 of the point difference satisfies a quadratic
        f01 over the base field.  Roots of f01 exist iff P and Q are on the same
        side; if they are on opposite sides the roots lie in the extension field.
        """
        self._check_same_surface(other, "other")
        K  = self.surface
        BB = K.biquadratic_matrix(self, other)
        R  = K.field['X']
        X  = R.gen()
        # f_ij(zi, X) = BB[j][j]*X^2 - 2*BB[i][j]*zi*X + BB[i][i]*zi^2
        # With Z0=1 affine, f01(Z0=1, Z1=X):
        f01 = BB[1][1]*X**2 - 2*BB[0][1]*X + BB[0][0]
        return len(f01.roots()) > 0

    def point_difference(self, other):
        """
        Compute (P - Q, P + Q) as a pair of SquaredKummerPoints, using the
        classical biquadratics on the intermediate Kummer surface (qDSA, Prop. 4).

        Raises ValueError if self and other lie on different Jacobians (the
        required roots then do not exist over the base field).
        """
        if not self.same_jacobian(other):
            raise ValueError(
                "lie on different Jacobians (Jacobian vs twist)."
            )
        
        self._check_same_surface(other, "other")
        K    = self.surface
        Kint = K.intermediate_kummer()
        Pint = K._hadamard(self.coords)
        Qint = K._hadamard(other.coords)
        BB   = K._biquadratic_matrix(Kint, Pint, Qint)

        # Work in the polynomial ring F[Z0, Z1, Z2, Z3].
        # The biquadratic form:
        #   F_ij = BB[i][i]*Z[j]^2 - 2*BB[i][j]*Z[i]*Z[j] + BB[j][j]*Z[i]^2
        R4 = K.field['Z0, Z1, Z2, Z3']
        Z  = list(R4.gens())
        F  = {(i, j): BB[i][i]*Z[j]**2 - 2*BB[i][j]*Z[i]*Z[j] + BB[j][j]*Z[i]**2
              for i in range(4) for j in range(i + 1, 4)}

        one = K.field(1)

        # Set Z0 = 1 (affine).  Z1 is a root of F_01 |_{Z0=1}.
        f01 = F[0, 1].subs({Z[0]: one}).univariate_polynomial()
        roots = [r for r, m in f01.roots() for _ in range(m)]

        if not roots:
            raise ValueError(
                "f01 has no roots over the base field: self and other likely "
                "lie on different Jacobians (Jacobian vs twist)."
            )

        diffs = []
        for z1 in roots:
            # Z2: common root of F_02 |_{Z0=1} and F_12 |_{Z1=z1}
            f02 = F[0, 2].subs({Z[0]: one}).univariate_polynomial()
            f12 = F[1, 2].subs({Z[1]: z1}).univariate_polynomial()
            g2  = f02.gcd(f12)
            z2  = -g2[0] / g2[1]

            # Z3: common root of F_03 |_{Z0=1} and F_13 |_{Z1=z1}
            f03 = F[0, 3].subs({Z[0]: one}).univariate_polynomial()
            f13 = F[1, 3].subs({Z[1]: z1}).univariate_polynomial()
            g3  = f03.gcd(f13)
            z3  = -g3[0] / g3[1]

            # Point (Z0=1 : Z1=z1 : Z2=z2 : Z3=z3) on intermediate surface;
            # map back to the squared Kummer surface via Hadamard.
            Dint = (one, z1, z2, z3)
            diffs.append(SquaredKummerPoint(K, K._hadamard(Dint)))

        S = diffs[0]
        D = diffs[1]
        
        assert self.xADD(other, S) == D
        assert self.xADD(other, D) == S
        assert S.xADD(D, 2*other) == 2*self
        
        return diffs[0], diffs[1]

    # ------------------------------------------------------------------
    # Doubling and differential addition
    # ------------------------------------------------------------------

    def hadamard(self):
        return SquaredKummerPoint(self.surface, self.surface._hadamard(self.coords))
    
    def sqr(self):
        return SquaredKummerPoint(self.surface, self.surface._k_sqr(self.coords))

    def __neg__(self):
        """this is of course ridiculous"""
        return self

    def xDBL(self):
        """Return 2 * self."""
        K = self.surface
        P = K._hadamard(self.coords)
        P = K._k_sqr(P)
        P = K._k_mul(P, K.K1)
        P = K._hadamard(P)
        P = K._k_sqr(P)
        P = K._k_mul(P, K.K0)
        return SquaredKummerPoint(K, P)

    def xDBLADD(self, Q, PQ):
        """
        Compute (2*self, self + Q) simultaneously.

        Q and PQ must be SquaredKummerPoints on the same surface.
        PQ must represent self - Q (equivalently self + Q on the Kummer).
        Returns a pair (doubled, added).
        """
        self._check_same_surface(Q, "Q")
        self._check_same_surface(PQ, "PQ")
        
        assert prod(PQ.coords) != 0

        K = self.surface
        HP = K._hadamard(self.coords)
        HQ = K._hadamard(Q.coords)

        U   = K._k_mul(HP, K.K1)          # K1 * H(P)
        dbl = K._k_mul(HP, U)             # K1 * H(P)^2
        add = K._k_mul(HQ, U)             # K1 * H(P) * H(Q)

        dbl = K._hadamard(dbl)
        add = K._hadamard(add)
        dbl = K._k_sqr(dbl)
        add = K._k_sqr(add)
        dbl = K._k_mul(dbl, K.K0)

        pq_inv = tuple(1/x for x in PQ.coords)
        add    = K._k_mul(add, pq_inv)

        return SquaredKummerPoint(K, dbl), SquaredKummerPoint(K, add)

    def xADD(self, Q, PQ):
        """
        Return self + Q given the difference PQ = self - Q.
        (Less efficient than xDBLADD when both 2*self and self+Q are needed.)
        """
        self._check_same_surface(Q, "Q")
        self._check_same_surface(PQ, "PQ")
        _, result = self.xDBLADD(Q, PQ)
        return result

    # ------------------------------------------------------------------
    # Scalar multiplication (Montgomery ladder)
    # ------------------------------------------------------------------

    def ladder(self, k):
        """Return k * self using a Montgomery-like ladder. k must be a positive integer."""
        if isinstance(k, bool):
            raise TypeError(f"Scalar k must be an integer, got {type(k)}")
        try:
            k = _int(k)
        except TypeError:
            raise TypeError(f"Scalar k must be an integer, got {type(k)}")
        if k <= 0:
            raise ValueError(f"Scalar k must be positive, got {k}")
        
        # i think this bug fix is enough
        if prod(self.coords) == 0 and k == 2:
            return self.surface.zero()

        K   = self.surface
        Pm  = SquaredKummerPoint(K, self.coords)
        Pp  = self.xDBL()
        for bit in bin(k)[3:]:   # skip leading '0b1'
            if bit == '0':
                Pm, Pp = Pm.xDBLADD(Pp, self)
            else:
                Pp, Pm = Pp.xDBLADD(Pm, self)
        return Pm

    def three_point_ladder(self, n, R, D):
        """
            from SIKE (for ell curves) and RotK (for Kummer surfaces)
            returns [n]*P + R
            using D = P - R
        """
        
        P0 = self
        P1 = R
        P2 = D

        for i in range(n.bit_length()):
            bit = (n >> i) & 1
            if bit == 1:
                P0, P1 = P0.xDBLADD(P1, P2)
            else:
                P0, P2 = P0.xDBLADD(P2, P1)
	
        return P1

    def __rmul__(self, k):
        """Support k * P syntax."""
        return self.ladder(k)

    def __mul__(self, k):
        """Support P * k syntax."""
        return self.ladder(k)

    # ------------------------------------------------------------------
    # Normalization and equality
    # ------------------------------------------------------------------

    def normalize(self):
        """Return a representative with the first non-zero coordinate set to 1."""
        for x in self.coords:
            if x != 0:
                ni = x
                break
        else:
            raise ValueError("Cannot normalize the all-zero tuple")
        return SquaredKummerPoint(
            self.surface,
            tuple(x / ni for x in self.coords),
            addition_matrix=self._addition_matrix,
        )

    def __eq__(self, other):
        if not isinstance(other, SquaredKummerPoint):
            return NotImplemented
        if other.surface is not self.surface:
            return False
        K = self.surface
        # Projective equality: P == Q iff k_mul(P, (Q[0],)*4) == k_mul(Q, (P[0],)*4)
        p0 = self.coords[0]
        q0 = other.coords[0]
        return K._k_mul(self.coords, (q0,)*4) == K._k_mul(other.coords, (p0,)*4)

    def __hash__(self):
        return hash(self.normalize().coords)

    def order(self):
        
        qo = self.surface.jacobian_order
        qt = self.surface.twist_order
        
        if (qo*self).is_zero():
            goal = qo
        elif (qt*self).is_zero():
            goal = qt            
        else:
            raise ValueError("Unclear situation in terms of order")

        res = 1
        for ell in factor(goal):
            Q = (goal // ell[0]**ell[1])*self
            while not Q.is_zero():
                res *= ell[0]
                Q *= ell[0]
                
        return res
    
    def order_factored(self):
        return factor(self.order())

    # ------------------------------------------------------------------
    # Membership tests
    # ------------------------------------------------------------------

    def is_zero(self):
        """Return True if this point equals the null point (identity)."""
        return self == self.surface.zero()

    def on_kummer(self):
        """Return True if this point satisfies the Kummer surface equation."""
        K = self.surface
        E, EE, F, G, H = K.curve
        X, Y, Z, T = self.coords
        LHS = EE * X*Y*Z*T
        mid = (X**2 + Y**2 + Z**2 + T**2
               - F*(X*T + Y*Z) - G*(X*Z + Y*T) - H*(X*Y + Z*T))
        RHS = mid**2
        return LHS == RHS

    def on_jacobian(self):
        """
        Return True if this point lies on the Jacobian (i.e. [jacobian_order]*self == 0).

        Requires the surface to have been constructed with a known jacobian_order.
        Raises NotImplementedError otherwise.
        """
        K = self.surface
        if K.jacobian_order is None:
            raise NotImplementedError(
                "jacobian_order is not set on this surface; "
                "pass jacobian_order= to SquaredKummerSurface() to enable this test."
            )
        return (K.jacobian_order * self).is_zero()

    def on_twist(self):
        """
        Return True if this point lies on the twist Jacobian (i.e. [twist_order]*self == 0).

        Requires the surface to have been constructed with a known twist_order.
        Raises NotImplementedError otherwise.
        """
        K = self.surface
        if K.twist_order is None:
            raise NotImplementedError(
                "twist_order is not set on this surface; "
                "pass twist_order= to SquaredKummerSurface() to enable this test."
            )
        return (K.twist_order * self).is_zero()
    
    def origin(self):
        if self.on_jacobian():
            return "Jacobian"
        elif self.on_twist():
            return "Twist"
        else:
            raise ValueError(
                "This should not happen: investigate closely."
            )
            
    def _invert(self):
        P = self.coords
        pi1 = P[2]*P[3]
        pi2 = pi1*P[0]
        pi1 = pi1*P[1]
        pi3 = P[0]*P[1]
        pi4 = pi3*P[2]
        pi3 = pi3*P[3]

        return [pi1,pi2,pi3,pi4]

    # ------------------------------------------------------------------
    # Pairing profile (cubical / Robert's algorithm)
    # ------------------------------------------------------------------

    def two_profile(self, basis=None):
        """
        Compute the 2-Tate pairing profile of this point.

        Returns a list of booleans: profile[i] is True iff
        t_2(basis[i], self) is a square in the base field.

        basis:  list of SquaredKummerPoints representing a basis of K[2].
                If None, the canonical basis from the surface is used
                (requires Rosenhain invariants on the surface).
        TODO: clean up, merge with other general profile computation
        """
        if basis is None:
            basis = self.surface.two_torsion_basis()
        if not isinstance(basis, (list, tuple)) or len(basis) == 0:
            raise ValueError("basis must be a non-empty list of SquaredKummerPoints")
        for T in basis:
            self._check_same_surface(T, "basis element")
        return [self._cubical_pairing_degree_2(T).is_square() for T in basis]

    def _cubical_pairing_degree_2(self, T):
        """
        Compute the reduced 2-Tate pairing t_2(T, self) via Robert's algorithm
        (Algorithm 5.2 / Appendix 3 of the paper).

        T must be a SquaredKummerPoint with _addition_matrix set.
        Returns a field element whose square-class encodes the pairing value.
        """
        K = self.surface
        M = T._addition_matrix
        if M is None:
            raise ValueError(
                "2-torsion basis points must have an associated addition matrix. "
                "Use surface.two_torsion_basis() to obtain them."
            )

        Q  = self.coords
        L  = K._matrix_mult(M, K._zero)       # L_ij = M * 0
        LQ = K._matrix_mult(M, Q)             # L_ij + Q = M * Q

        # Find a non-zero coordinate index of L
        ni = next(
            (i for i, x in enumerate(L) if x != 0),
            None,
        )
        if ni is None:
            raise ArithmeticError("2-torsion point L is the zero tuple; bad addition matrix")

        LT  = K._matrix_mult(M, L)            # M * L_ij
        LQT = K._matrix_mult(M, LQ)           # M * (L_ij + Q)

        lambda_L = LT[ni]  / L[ni]
        lambda_Q = LQT[ni] / LQ[ni]

        return lambda_Q / lambda_L
    
    # ------------------------------------------------------------------
    # HIGHER DEGREE PAIRINGS
    # ------------------------------------------------------------------

    def unreduced_tate_pairing(self, R, ell, difference=None):
        # essentially computes the square of the Tate pairing
        # TODO: fix linearity
        assert (ell*self).is_zero()
        P = self
        
        q = self.surface.field.cardinality()
        if (q - 1) % ell != 0:
            k = self.surface.embedding_degree(ell)
            K_ext = self.surface.base_change(k)    
            P = K_ext(self.coords) 
            R = K_ext(R.coords)
        
        if difference == None:
            S, D = P.point_difference(R)
        else:
            D = difference
        
        # P = P.normalize()
        # R = R.normalize()
        # D = D.normalize()
        
        P_ell_R = P.three_point_ladder(ell , R, D)
        assert P_ell_R == R
        
        lam1 = P_ell_R[0]/R[0]
        lam0 = (ell*P)[0]/P.surface._zero[0]
        
        return lam1/lam0

    def tate_pairing(self, R, ell, difference=None):
        # TODO: generalize to any field, not just Fp2
        assert ell % 2 == 1
        assert (ell*self).is_zero()
        q = self.surface.field.cardinality()
        k = self.surface.embedding_degree(ell)
        
        # TODO: final reduction should use frobenius
        return self.unreduced_tate_pairing(R, ell, difference=difference)**((q**k-1) // ell)
    
    def profile(self, ell):
        # TODO: generalize to even degrees too
        # when computed from Jacobian, can include differences
        assert ell % 2 == 1
        return TateProfile(tuple([ P.tate_pairing(self, ell) for P in self.surface.torsion_basis(ell) ]))

    def weil_pairing(self, other, ell, difference=None):
        assert (ell*self).is_zero()
        assert (ell*other).is_zero()
        res = self.unreduced_tate_pairing(other, ell, difference=difference) / other.unreduced_tate_pairing(self, ell, difference=difference)
        assert res**ell == 1
        return res
    
    # ------------------------------------------------------------------
    # Representation
    # ------------------------------------------------------------------

    def __repr__(self):
        x1, x2, x3, x4 = self.coords
        return f"({x1} : {x2} : {x3} : {x4})"

    def __getitem__(self, key):
        return self.coords[key]