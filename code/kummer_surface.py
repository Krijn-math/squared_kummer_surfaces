"""
SquaredKummerSurface: a Kummer surface defined by its null point.

A Kummer surface K is defined by four projective coordinates (a : b : c : d)
of its null point (the image of the identity of the Jacobian). All arithmetic
constants are derived from these four values.

Optionally, Rosenhain invariants (lam, mu, nu) can be provided to enable
computation of the 2-torsion addition matrices and the pairing profile.
"""

from sage.all import *


def hadamard(P):
    X, Y, Z, T = P
    return (X+Y+Z+T, X+Y-Z-T, X-Y+Z-T, X-Y-Z+T)


def k_sqr(P):
    return tuple(x**2 for x in P)


def k_mul(P, Q):
    return tuple(p * q for p, q in zip(P, Q))


class SquaredKummerSurface:

    def __init__(self, zero, jacobian_order=None, twist_order=None):
        """
        Initialize a Kummer surface from its null point.

        zero:            4-tuple (a, b, c, d) of elements in a Sage finite field.
        jacobian_order:  optional integer order of the Jacobian above this Kummer surface.
                         Required for SquaredKummerPoint.on_jacobian().
        twist_order:     optional integer order of the twist Jacobian above this Kummer surface.
                         Required for SquaredKummerPoint.on_twist().
        """
        
        ## when defined by null point
        if len(zero) == 4:
            self._zero = tuple(zero)
        
        # when defined by rosenhain invariants
        elif len(zero) == 3:
            ros = zero
            self.rosenhain = tuple(ros)
            self._zero = self.rosenhain_to_zero()
            
        self.field = zero[0].parent()
        self._compute_constants()
                    
        try:
            self.rosenhain
        except:
            self.rosenhain = self.rosenhain_invariants()    

        self.jacobian_order = jacobian_order
        self.twist_order = twist_order
        self.Mij = None
        self._torsion_basis_cache = {}
        self._compute_addition_matrices()
        
    def mumford_to_squared_theta(self, point):
        # assuming that the mumford representation
        # is given for a curve in rosenhain form
        
        λ, μ, ν = point.scheme()._ros
        u0, u1, _ = point[0].coefficients()
        v0, _ = point[1].coefficients()

        # K = SquaredKummerSurface(ros)
        # a, b, c, d = K.zero()
        a, b, c, d = self.zero()

        X1 = a * (u0 * (1*μ - u0) * (λ + ν + u1) - v0**2)
        X2 = b * (u0 * (λ*ν - u0) * (1 + μ + u1) - v0**2)
        X3 = c * (u0 * (1*ν - u0) * (λ + μ + u1) - v0**2)
        X4 = d * (u0 * (λ*μ - u0) * (1 + ν + u1) - v0**2)
        
        return (X1, X2, X3, X4)

    def __call__(self, point):
        from kummer_point import SquaredKummerPoint
        
        if isinstance(point, SquaredKummerPoint):
            assert point.on_kummer()
            assert point.surface == self
            return point
        
        if isinstance(point, tuple):
            return SquaredKummerPoint(self, point)
        
        if isinstance(point[0], Polynomial):
            return SquaredKummerPoint(self, self.mumford_to_squared_theta(point))

        raise NotImplementedError("Unclear what input is")
    
    def __eq__(self, other):
        p0 = self.zero()[0]
        q0 = other.zero()[0]
        return self._k_mul(self._zero, (q0,)*4) == other._k_mul(other._zero, (p0,)*4)
    
    def zero(self):
        from kummer_point import SquaredKummerPoint
        return SquaredKummerPoint(self, self._zero, self.Mij)
        
    def base_change(self, degree):
        Fk = self.field.extension(degree)
        new_zero = tuple([ Fk(z) for z in self._zero ])
        return SquaredKummerSurface(new_zero, jacobian_order=self.jacobian_order, twist_order=self.twist_order)
    
    def embedding_degree(self, ell):
        """
            returns embedding degree k as smallest value such that ell divides q^k -1
        """
        q = self.field.cardinality()
        k = 2
        while (q**k - 1) % ell != 0:
            k += 1
        
        return k

    # ------------------------------------------------------------------
    # Internal arithmetic helpers (operate on plain tuples)
    # ------------------------------------------------------------------

    def _hadamard(self, P):
        return hadamard(P)

    def _k_sqr(self, P):
        return k_sqr(P)

    def _k_mul(self, P, Q):
        return k_mul(P, Q)

    def _matrix_mult(self, M, P):
        """Compute M * P where M is a 4x4 list-of-lists and P is a 4-tuple."""
        return tuple(sum(M[j][i] * P[i] for i in range(4)) for j in range(4))

    # ------------------------------------------------------------------
    # Derived constants
    # ------------------------------------------------------------------

    def _compute_constants(self):
        a, b, c, d = self._zero

        # Dual null point via Hadamard transform
        self.dual_zero = hadamard(self._zero)
        A, B, C, D = self.dual_zero

        # Componentwise inverse of null point
        self.zero_inv = tuple(1/x for x in self._zero)

        # Doubling constants:
        #   K0[i] = zero[0] / zero[i]   (applied after second square)
        #   K1[i] = dual_zero[0] / dual_zero[i]  (applied after first square)
        self.K0 = (self.field(1), a/b, a/c, a/d)
        self.K1 = (self.field(1), A/B, A/C, A/D)

        # Kummer surface equation constants
        # (See Cassels-Flynn or the paper; computed from the null point.)
        zero2 = k_sqr(self._zero)          # (a^2, b^2, c^2, d^2)
        a2, b2, c2, d2 = a, b, c, d
        a4, b4, c4, d4 = zero2

        E_top = A * B * C * D
        E_bot = (a2*d2 - b2*c2) * (a2*c2 - b2*d2) * (a2*b2 - c2*d2)
        E  = E_top / E_bot
        EE = E**2 * 4 * a2 * b2 * c2 * d2

        F = (a4 - b4 - c4 + d4) / (a2*d2 - b2*c2)
        G = (a4 - b4 + c4 - d4) / (a2*c2 - b2*d2)
        H = (a4 + b4 - c4 - d4) / (a2*b2 - c2*d2)

        self.curve = (E, EE, F, G, H)

    # ------------------------------------------------------------------
    # Rosenhain invariants
    # ------------------------------------------------------------------

    def rosenhain_to_zero(self):
        ros = self.rosenhain
        λ, μ, ν = ros
        d = 1
        c = (λ*μ/ν).sqrt()
        b = ( (μ*(μ-1)*(λ-ν)) / (ν*(ν-1)*(λ - μ)) ).sqrt()
        a = b*c*ν/μ
        return (a, b, c, d)

    def rosenhain_invariants(self, bit=0):
        """
        Compute the Rosenhain invariants (lam, mu, nu) for this Kummer surface.

        If Rosenhain invariants were already provided at construction, they are
        returned directly (bit is ignored).

        Otherwise they are derived from the null point via one square root.
        Given null point (a,b,c,d) with dual (A,B,C,D) = hadamard(a,b,c,d):

            alpha^2 = C*D / (A*B)
            e/f     = (1+alpha) / (1-alpha)
            lam     = a*c / (b*d)
            mu      = (c/d) * (e/f)
            nu      = (a/b) * (e/f)

        The ratio e/f (and hence all three invariants) depends only on the null
        point.  bit=0 uses alpha = +sqrt(C*D/(A*B)); bit=1 uses -alpha,
        selecting the conjugate Rosenhain triple (Jacobian vs twist).

        Returns a 3-tuple (lam, mu, nu) of field elements.
        """
        try:
            return self.rosenhain
        except:
            a, b, c, d = self._zero
            A, B, C, D = self.dual_zero

            alpha = (C * D / (A * B)).sqrt()
            if bit:
                alpha = -alpha

            ratio = (1 + alpha) / (1 - alpha)

            lam = a * c / (b * d)
            mu  = c * ratio / d
            nu  = a * ratio / b

            return (lam, mu, nu)

    # ------------------------------------------------------------------
    # 2-torsion addition matrices
    # ------------------------------------------------------------------

    def _compute_addition_matrices(self):
        """
        Compute the 4x4 addition matrices Mij[i][j] (0 <= i < j <= 5) from
        the Rosenhain invariants and null point.

        Each matrix Mij[i][j] represents translation by the 2-torsion point
        L_{ij} on the Kummer surface, i.e. it maps P -> P + L_{ij}.

        The formula follows matrices_precompute_script.py, derived from the
        Rosenhain model of the underlying genus-2 curve.
        """
        F = self.field
        lam, mu, nu = self.rosenhain
        a, b, c, d = self._zero

        # Weierstrass points of the Rosenhain model y^2 = x(x-1)(x-lam)(x-mu)(x-nu):
        #   w[0] = 0, w[1] = 0, w[2] = 0, w[3] = 1, w[4] = lam, w[5] = mu, w[6] = nu
        # (indices 0..6; only indices 2..6 appear in the formulas)
        w = [F(0), F(0), F(0), F(1), lam, mu, nu]

        # Null point coordinates as mus[1..4]
        m = [F(0), a, b, c, d]

        Z = [[F(0)]*4 for _ in range(4)]
        Mij = [[None]*6 for _ in range(6)]

        # Helper: build a 4-tuple row from four values
        def row(*vals):
            return [F(v) for v in vals]

        # --- Permutation matrices (indices whose formulas degenerate) ---
        Mij[0][1] = [row(0,1,0,0), row(1,0,0,0), row(0,0,0,1), row(0,0,1,0)]
        Mij[2][3] = [row(0,0,0,1), row(0,0,1,0), row(0,1,0,0), row(1,0,0,0)]
        Mij[4][5] = [row(0,0,1,0), row(0,0,0,1), row(1,0,0,0), row(0,1,0,0)]

        # --- General matrices computed from Rosenhain formulas ---
        # Shorthand: r(i,j) = wp[i] - wp[j]
        def r(i, j):
            return w[i] - w[j]

        # M[0][2]
        alpha = r(2,5)/r(3,2) * m[4]/m[3]
        beta  = r(2,5)/r(2,3) * m[4]/m[3]
        gamma = r(4,5)/r(5,3) * m[2]/m[3]
        delta = r(4,5)/r(3,5) * m[2]/m[3]
        eps   = beta / (r(3,5)/r(4,5) * m[3]/m[2])
        zeta  = beta / (r(3,5)/r(5,4) * m[3]/m[2])
        Mij[0][2] = [
            row(1,    alpha, eps,   gamma),
            row(beta, -1,    delta, zeta ),
            row(eps,  gamma, 1,     alpha),
            row(delta,zeta,  beta,  -1   ),
        ]

        # M[0][3]
        alpha = r(2,5)/r(3,2) * m[4]/m[3]
        beta  = r(2,5)/r(2,3) * m[4]/m[3]
        gamma = r(3,5)/r(5,4) * m[3]/m[2]
        delta = r(3,5)/r(4,5) * m[3]/m[2]
        eps   = beta / (r(4,5)/r(3,5) * m[2]/m[3])
        zeta  = beta / (r(4,5)/r(5,3) * m[2]/m[3])
        Mij[0][3] = [
            row(1,    alpha, eps,   gamma),
            row(beta, -1,    delta, zeta ),
            row(eps,  gamma, 1,     alpha),
            row(delta,zeta,  beta,  -1   ),
        ]

        # M[0][4]
        alpha = r(2,3)/r(5,2) * m[3]/m[4]
        beta  = r(2,3)/r(2,5) * m[3]/m[4]
        gamma = r(3,6)/r(5,3) * m[2]/m[4]
        delta = r(3,6)/r(3,5) * m[2]/m[4]
        eps   = beta / (r(3,5)/r(3,6) * m[4]/m[2])
        zeta  = beta / (r(3,5)/r(6,3) * m[4]/m[2])
        Mij[0][4] = [
            row(1,    alpha, gamma, eps  ),
            row(beta, -1,    zeta,  delta),
            row(delta,zeta,  -1,    beta ),
            row(eps,  gamma, alpha, 1    ),
        ]

        # M[0][5]
        alpha = r(2,3)/r(5,2) * m[3]/m[4]
        beta  = r(2,3)/r(2,5) * m[3]/m[4]
        gamma = r(3,5)/r(6,3) * m[4]/m[2]
        delta = r(3,5)/r(3,6) * m[4]/m[2]
        eps   = beta / (r(3,6)/r(3,5) * m[2]/m[4])
        zeta  = beta / (r(3,6)/r(5,3) * m[2]/m[4])
        Mij[0][5] = [
            row(1,    alpha, gamma, eps  ),
            row(beta, -1,    zeta,  delta),
            row(delta,zeta,  -1,    beta ),
            row(eps,  gamma, alpha, 1    ),
        ]

        # M[1][2]
        alpha = r(2,3)/r(5,2) * m[3]/m[4]
        beta  = r(2,3)/r(2,5) * m[3]/m[4]
        gamma = r(4,5)/r(5,3) * m[2]/m[3]
        delta = r(4,5)/r(3,5) * m[2]/m[3]
        eps   = beta / (r(3,5)/r(4,5) * m[3]/m[2])
        zeta  = beta / (r(3,5)/r(5,4) * m[3]/m[2])
        Mij[1][2] = [
            row(1,    alpha, eps,   gamma),
            row(beta, -1,    delta, zeta ),
            row(eps,  gamma, 1,     alpha),
            row(delta,zeta,  beta,  -1   ),
        ]

        # M[1][3]
        alpha = r(2,3)/r(5,2) * m[3]/m[4]
        beta  = r(2,3)/r(2,5) * m[3]/m[4]
        gamma = r(3,5)/r(5,4) * m[3]/m[2]
        delta = r(3,5)/r(4,5) * m[3]/m[2]
        eps   = beta / (r(4,5)/r(3,5) * m[2]/m[3])
        zeta  = beta / (r(4,5)/r(5,3) * m[2]/m[3])
        Mij[1][3] = [
            row(1,    alpha, eps,   gamma),
            row(beta, -1,    delta, zeta ),
            row(eps,  gamma, 1,     alpha),
            row(delta,zeta,  beta,  -1   ),
        ]

        # M[1][4]
        alpha = r(2,5)/r(3,2) * m[4]/m[3]
        beta  = r(2,5)/r(2,3) * m[4]/m[3]
        gamma = r(3,6)/r(5,3) * m[2]/m[4]
        delta = r(3,6)/r(3,5) * m[2]/m[4]
        eps   = beta / (r(3,5)/r(3,6) * m[4]/m[2])
        zeta  = beta / (r(3,5)/r(6,3) * m[4]/m[2])
        Mij[1][4] = [
            row(1,    alpha, gamma, eps  ),
            row(beta, -1,    zeta,  delta),
            row(delta,zeta,  -1,    beta ),
            row(eps,  gamma, alpha, 1    ),
        ]

        # M[1][5]
        alpha = r(2,5)/r(3,2) * m[4]/m[3]
        beta  = r(2,5)/r(2,3) * m[4]/m[3]
        gamma = r(3,5)/r(6,3) * m[4]/m[2]
        delta = r(3,5)/r(3,6) * m[4]/m[2]
        eps   = beta / (r(3,6)/r(3,5) * m[2]/m[4])
        zeta  = beta / (r(3,6)/r(5,3) * m[2]/m[4])
        Mij[1][5] = [
            row(1,    alpha, gamma, eps  ),
            row(beta, -1,    zeta,  delta),
            row(delta,zeta,  -1,    beta ),
            row(eps,  gamma, alpha, 1    ),
        ]

        # M[2][4]
        alpha = r(3,5)/r(4,6) * m[1]/m[2]
        beta  = r(3,5)/r(6,3) * m[4]/m[2]
        gamma = r(3,5)/r(5,4) * m[3]/m[2]
        delta = r(3,5)/r(3,6) * m[4]/m[2]
        eps   = r(3,5)/r(4,5) * m[3]/m[2]
        zeta  = r(3,5)/r(6,4) * m[1]/m[2]
        Mij[2][4] = [
            row(1,     alpha, beta,  gamma),
            row(alpha, 1,     gamma, beta ),
            row(delta, eps,   -1,    zeta ),
            row(eps,   delta, zeta,  -1   ),
        ]

        # M[2][5]
        alpha = r(3,6)/r(4,5) * m[3]/m[4]
        beta  = r(3,6)/r(5,3) * m[2]/m[4]
        gamma = r(3,5)/r(5,4) * m[3]/m[2]
        delta = r(3,6)/r(3,5) * m[2]/m[4]
        eps   = r(3,5)/r(4,5) * m[3]/m[2]
        zeta  = r(3,6)/r(5,4) * m[3]/m[4]
        Mij[2][5] = [
            row(1,     alpha, beta,  gamma),
            row(alpha, 1,     gamma, beta ),
            row(delta, eps,   -1,    zeta ),
            row(eps,   delta, zeta,  -1   ),
        ]

        # M[3][4]
        alpha = r(4,5)/r(3,6) * m[4]/m[3]
        beta  = r(3,5)/r(6,3) * m[4]/m[2]
        gamma = r(4,5)/r(5,3) * m[2]/m[3]
        delta = r(3,5)/r(3,6) * m[4]/m[2]
        eps   = r(4,5)/r(3,5) * m[2]/m[3]
        zeta  = r(4,5)/r(6,3) * m[4]/m[3]
        Mij[3][4] = [
            row(1,     alpha, beta,  gamma),
            row(alpha, 1,     gamma, beta ),
            row(delta, eps,   -1,    zeta ),
            row(eps,   delta, zeta,  -1   ),
        ]

        # M[3][5]
        alpha = r(4,6)/r(3,5) * m[2]/m[1]
        beta  = r(3,6)/r(5,3) * m[2]/m[4]
        gamma = r(4,5)/r(5,3) * m[2]/m[3]
        delta = r(3,6)/r(3,5) * m[2]/m[4]  # wait, need to re-check
        eps   = r(4,5)/r(3,5) * m[2]/m[3]
        zeta  = r(4,6)/r(5,3) * m[2]/m[1]
        Mij[3][5] = [
            row(1,     alpha, beta,  gamma),
            row(alpha, 1,     gamma, beta ),
            row(delta, eps,   -1,    zeta ),
            row(eps,   delta, zeta,  -1   ),
        ]

        self.Mij = Mij

    # ------------------------------------------------------------------
    # Point constructors
    # ------------------------------------------------------------------

    def null_point(self):
        """Return the null point (identity) as a SquaredKummerPoint."""
        from kummer_point import SquaredKummerPoint
        return SquaredKummerPoint(self, self._zero)

    def point(self, coords, validate=True):
        """
        Create a SquaredKummerPoint on this surface.

        validate:  if True (default), raise ValueError if the coordinates do not
                   satisfy the Kummer surface equation.
        """
        from kummer_point import SquaredKummerPoint
        pt = SquaredKummerPoint(self, tuple(coords))
        if validate and not pt.on_kummer():
            raise ValueError(
                f"Coordinates {tuple(coords)} do not lie on the Kummer surface."
            )
        return pt

    def random_point(self, point=None):
        """
        Return a random point on this Kummer surface.

        Picks random values X, Y, Z from the field and solves the degree-4
        Kummer equation for T.  Retries until a solution exists.
        """
        from kummer_point import SquaredKummerPoint
        F = self.field
        _, EE, Fv, Gv, Hv = self.curve
        R = F['T']
        T_var = R.gen()

        while True:
            X = F.random_element()
            Y = F.random_element()
            Z = F.random_element()
            if X == 0 or Y == 0 or Z == 0:
                continue

            c1 = -Fv*X - Gv*Y - Hv*Z
            c0 = X**2 + Y**2 + Z**2 - Fv*Y*Z - Gv*X*Z - Hv*X*Y

            mid = T_var**2 + c1*T_var + c0
            poly = mid**2 - EE*X*Y*Z*T_var

            roots = poly.roots()
            if roots:
                T = roots[0][0]
                res = SquaredKummerPoint(self, (X, Y, Z, T))
                if point == None:
                    return res
                elif point.same_jacobian(res):
                    return res

    def two_torsion_basis(self):
        """
        Return the canonical 2-torsion basis as a list of four SquaredKummerPoints.

        The basis consists of the 2-torsion points L_{ij} for indices
        (i,j) in {(1,2), (1,3), (1,4), (1,5)}, derived from the Rosenhain model.

        Requires Rosenhain invariants to have been provided at construction.
        """
        from kummer_point import SquaredKummerPoint
        if self.Mij is None:
            raise ValueError(
                "Rosenhain invariants are required to compute the 2-torsion basis."
            )
        basis_indices = [(1, 2), (1, 3), (1, 4), (1, 5)]
        basis = []
        for i, j in basis_indices:
            M = self.Mij[i][j]
            coords = self._matrix_mult(M, self._zero)
            pt = SquaredKummerPoint(self, coords, addition_matrix=M)
            basis.append(pt)
        return basis

    def slow_multiples(self, ell, point):
        # TODO: improve this for comparison for basis
        return [k*point for k in range(1,(ell + 1) // 2)]

    def torsion_basis(self, ell):
        """
            for now assumes ell divides p+1 and supersingular kummer surfaces
        """
        if ell in self._torsion_basis_cache:
            return self._torsion_basis_cache[ell]

        p = self.field.characteristic()
        assert (p+1) % ell == 0
            
        # first point on Jac
        P1 = self.point_of_order(ell)
        while True:
            # second different point
            P2 = self.point_of_order(ell)
            # TODO: improve check for multiples
            if P1.same_jacobian(P2) and P2 not in self.slow_multiples(ell, P1):
                break
            
        # then point on twist
        while True:
            Q1 = self.point_of_order(ell)
            if not P1.same_jacobian(Q1):
                break
            
        while True:
            Q2 = self.point_of_order(ell)
            if Q1.same_jacobian(Q2) and Q2 not in self.slow_multiples(ell, Q1):
                break
            
        basis = (P1, P2, Q1, Q2)
        self._torsion_basis_cache[ell] = basis
        return basis

    def __repr__(self):
        return f"Squared Kummer Surface defined by zero: {self.zero()}"
    
    def point_of_order(self, ell):
        if (self.jacobian_order % ell *self.twist_order % ell) != 0:
            return "Order seems infeasible in this base field"

        while True:
            P = self.random_point()
            if P.order() % ell == 0:
                break
            
        P = (P.order() // ell)*P
        assert not P.is_zero()
        return P
    
    def intermediate_kummer(self):
        """
        Compute constants for the intermediate Kummer surface (cached).

        Following qDSA (Renes-Smith), the intermediate surface uses
        mudual = (1/2) * hadamard(null_point) and a scaling constant C.

        Returns (mu, mudual, C) where mu is the null point, mudual its
        scaled dual, and C a scalar normalization constant.
        """
        if hasattr(self, '_kint'):
            return self._kint
        mu = self._zero
        inv2 = self.field(1) / 2
        mudual = tuple(inv2 * m for m in hadamard(mu))
        A, B, C, D = mudual
        C1 = A*B - C*D
        C2 = A*C - B*D
        C3 = A*D - B*C
        Cv = mu[0]*mu[1]*mu[2]*mu[3] * (A*B*C*D) * 8 / (C1 * C2 * C3)
        self._kint = (mu, mudual, Cv)
        return self._kint

    def biquadratic_matrix(self, P, Q):
        """
        Compute the 4x4 symmetric biquadratic matrix for SquaredKummerPoints P, Q.

        P and Q are mapped to the intermediate Kummer surface (Hadamard applied)
        before evaluation.  Returns a 4x4 list of field elements.
        """
        Kint = self.intermediate_kummer()
        Pint = self._hadamard(P.coords)
        Qint = self._hadamard(Q.coords)
        return self._biquadratic_matrix(Kint, Pint, Qint)

    def _biquadratic_matrix(self, Kint, Pint, Qint):
        """Internal: Pint, Qint are raw tuples already on the intermediate surface."""
        BB = [[None]*4 for _ in range(4)]
        for i in range(4):
            BB[i][i] = self._bii(Kint, i, Pint, Qint)
        for i in range(3):
            for j in range(i + 1, 4):
                v = self._bij(Kint, i, j, Pint, Qint)
                BB[i][j] = v
                BB[j][i] = v
        return BB

    def _bii(self, Kint, i, P, Q):
        """Diagonal (i,i) entry of the biquadratic matrix."""
        mu, mudual, _Cv = Kint
        mu01 = mu[0]*mu[1]
        mu23 = mu[2]*mu[3]
        eps    = (mu[1]*mu23, mu[0]*mu23, mu01*mu[3], mu01*mu[2])
        kappa  = hadamard(eps)
        
        def invert_constants(P):
            pi1 = P[2]*P[3]
            pi2 = pi1*P[0]
            pi1 = pi1*P[1]
            pi3 = P[0]*P[1]
            pi4 = pi3*P[2]
            pi3 = pi3*P[3]

            return [pi1,pi2,pi3,pi4]

        epsinv = invert_constants(mudual)
        Pe = k_mul(k_sqr(P), epsinv)
        Qe = k_mul(k_sqr(Q), epsinv)
        F1 = sum(k_mul(Pe, Qe))
        F2 = sum(k_mul(Pe, (Qe[1], Qe[0], Qe[3], Qe[2])))
        F3 = sum(k_mul(Pe, (Qe[2], Qe[3], Qe[0], Qe[1])))
        F4 = sum(k_mul(Pe, (Qe[3], Qe[2], Qe[1], Qe[0])))
        # kappa index permutation for each i
        perm = ((0,1,2,3), (1,0,3,2), (2,3,0,1), (3,2,1,0))
        ks = perm[i]
        return mudual[i] * (kappa[ks[0]]*F1 + kappa[ks[1]]*F2
                          + kappa[ks[2]]*F3 + kappa[ks[3]]*F4)

    def _bij(self, Kint, i, j, P, Q):
        """Off-diagonal (i,j) entry of the biquadratic matrix."""
        _mu, mudual, Cv = Kint
        k, l = sorted({0, 1, 2, 3} - {i, j})
        C1  = mudual[i]*mudual[j]
        C2  = mudual[i]*mudual[k] - mudual[j]*mudual[l]
        C3  = mudual[i]*mudual[l] - mudual[j]*mudual[k]
        Cij = C1 * C2 * C3
        Pij = P[i]*P[j];  Pkl = P[k]*P[l]
        Qij = Q[i]*Q[j];  Qkl = Q[k]*Q[l]
        C4      = mudual[i]*mudual[j] - mudual[k]*mudual[l]
        Bij_pre = (Pij - Pkl) * (Qij - Qkl) * mudual[k] * mudual[l]
        Bij_pre += C4 * Pkl * Qkl
        return Cv * Cij * Bij_pre