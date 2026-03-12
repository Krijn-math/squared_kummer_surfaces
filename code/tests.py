"""
Tests for SquaredKummerSurface and SquaredKummerPoint.

Run with:  sage -python -m pytest tests.py
       or: python -m pytest tests.py  (if sage is on the path)
"""

import unittest
from random import choice

from example import K, gen, qo, random_points, basis, F


class TestKummerEquation(unittest.TestCase):
    """Points that should satisfy the Kummer surface equation."""

    def test_null_point_on_kummer(self):
        self.assertTrue(K.null_point().on_kummer())

    def test_generator_on_kummer(self):
        self.assertTrue(gen.on_kummer())

    def test_random_points_on_kummer(self):
        for P in random_points[:10]:
            self.assertTrue(P.on_kummer(), f"Point not on Kummer: {P}")


class TestDoubling(unittest.TestCase):
    """Doubling the null point returns the null point."""

    def test_dbl_null_is_null(self):
        Z = K.null_point()
        self.assertEqual(Z.xDBL().normalize(), Z.normalize())

    def test_dbl_random(self):
        P = choice(random_points)
        self.assertTrue(P.xDBL().on_kummer())


class TestLadder(unittest.TestCase):
    """Scalar multiplication correctness."""

    def test_order_of_generator(self):
        self.assertTrue((qo * gen).is_zero(), "qo * gen should be zero")

    def test_double_via_ladder(self):
        # 2*P should equal P.xDBL()
        P = choice(random_points)
        self.assertEqual((2 * P).normalize(), P.xDBL().normalize())

    def test_rmul_and_mul_agree(self):
        P = choice(random_points)
        k = 17
        self.assertEqual((k * P).normalize(), (P * k).normalize())


class TestProfile(unittest.TestCase):
    """Pairing profile properties."""

    def _trivial_profile(self):
        return [True, True, True, True]

    def test_profile_of_doubled_point_is_trivial(self):
        # A doubled point is in the image of [2]: Jac -> Jac,
        # so its profile against any 2-torsion basis should be trivial.
        for P in random_points[:10]:
            prof = P.xDBL().profile(basis)
            self.assertEqual(prof, self._trivial_profile(),
                             f"Non-trivial profile for 2P: {P}")

    def test_profile_returns_list_of_booleans(self):
        P = choice(random_points)
        prof = P.profile(basis)
        self.assertIsInstance(prof, list)
        self.assertEqual(len(prof), 4)
        for v in prof:
            self.assertIsInstance(v, bool)

    def test_profile_default_basis(self):
        # profile() with no argument should give same result as explicit basis
        P = choice(random_points)
        self.assertEqual(P.profile(), P.profile(basis))


class TestTwoTorsionBasis(unittest.TestCase):
    """The canonical 2-torsion basis elements satisfy [2]*T == 0 and T != 0."""

    def test_basis_elements_have_order_two(self):
        for T in basis:
            self.assertFalse(T.is_zero(), "Basis element should not be zero")
            self.assertTrue(T.xDBL().is_zero(),
                            f"2*T should be zero, got {T.xDBL()}")

    def test_basis_elements_on_kummer(self):
        for T in basis:
            self.assertTrue(T.on_kummer())


class TestRandomPoint(unittest.TestCase):
    """Random point generation on the Kummer surface."""

    def test_random_point_on_kummer(self):
        for _ in range(5):
            P = K.random_point()
            self.assertTrue(P.on_kummer(), f"Random point not on Kummer: {P}")

    def test_random_point_is_kummer_point(self):
        from kummer_point import SquaredKummerPoint
        P = K.random_point()
        self.assertIsInstance(P, SquaredKummerPoint)
        self.assertIs(P.surface, K)

    def test_random_points_are_distinct(self):
        pts = [K.random_point() for _ in range(3)]
        # Normalise and compare; collision is astronomically unlikely
        normed = [p.normalize().coords for p in pts]
        self.assertEqual(len(set(normed)), 3)


class TestSubgroupMembership(unittest.TestCase):
    """on_jacobian() and on_twist() membership tests."""

    def test_generator_on_jacobian(self):
        self.assertTrue(gen.on_jacobian())

    def test_generator_not_on_twist(self):
        self.assertFalse(gen.on_twist())

    def test_on_jacobian_raises_without_order(self):
        from kummer_surface import SquaredKummerSurface
        K2 = SquaredKummerSurface(K.zero)
        P = K2.null_point()
        with self.assertRaises(NotImplementedError):
            P.on_jacobian()

    def test_on_twist_raises_without_order(self):
        from kummer_surface import SquaredKummerSurface
        K2 = SquaredKummerSurface(K.zero)
        P = K2.null_point()
        with self.assertRaises(NotImplementedError):
            P.on_twist()


class TestPointValidation(unittest.TestCase):
    """surface.point() validates coordinates against the Kummer equation."""

    def test_valid_coords_accepted(self):
        # null point coordinates are always on the surface
        P = K.point(K.zero)
        self.assertTrue(P.on_kummer())

    def test_invalid_coords_rejected(self):
        bad = (K.field(1), K.field(2), K.field(3), K.field(4))
        with self.assertRaises(ValueError):
            K.point(bad)

    def test_validate_false_skips_check(self):
        bad = (K.field(1), K.field(2), K.field(3), K.field(4))
        # Should not raise even for invalid coords
        P = K.point(bad, validate=False)
        self.assertFalse(P.on_kummer())


class TestHygiene(unittest.TestCase):
    """Input validation."""

    def test_rmul_rejects_float(self):
        P = choice(random_points)
        with self.assertRaises(TypeError):
            _ = 1.0 * P

    def test_rmul_rejects_bool(self):
        P = choice(random_points)
        with self.assertRaises(TypeError):
            _ = True * P

    def test_ladder_rejects_zero(self):
        P = choice(random_points)
        with self.assertRaises(ValueError):
            P.ladder(0)

    def test_ladder_rejects_negative(self):
        P = choice(random_points)
        with self.assertRaises(ValueError):
            P.ladder(-1)

    def test_xdbladd_rejects_wrong_surface(self):
        from kummer_surface import SquaredKummerSurface
        other_zero = (F(11), F(-22), F(-19), F(-4))
        K2 = SquaredKummerSurface(other_zero)
        Q = K2.null_point()
        P = choice(random_points)
        with self.assertRaises(ValueError):
            P.xDBLADD(Q, Q)

    def test_profile_rejects_empty_basis(self):
        P = choice(random_points)
        with self.assertRaises(ValueError):
            P.profile([])


if __name__ == "__main__":
    unittest.main()
