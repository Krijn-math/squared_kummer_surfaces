"""
Example: the specific Kummer surface over GF(2^127 - 1).

This instantiates the SquaredKummerSurface and SquaredKummerPoint classes
for the Kummer surface used in the pairing research.  All hardcoded constants
(null point, Rosenhain invariants, curve orders, generator, random test points)
live here so the library classes remain parameter-free.
"""

from sage.all import GF
from kummer_surface import SquaredKummerSurface
from kummer_point import SquaredKummerPoint
from points import randompoints as _raw_points

# ------------------------------------------------------------------
# Base field
# ------------------------------------------------------------------
p = (1 << 127) - 1
F = GF(p)

# ------------------------------------------------------------------
# Kummer surface: null point and Rosenhain invariants
# ------------------------------------------------------------------
_zero = (F(11), F(-22), F(-19), F(-3))

_lam = F(28356863910078205288614550619314017618)
_mu  = F(154040945529144206406682019582013187910)
_nu  = F(113206060534360680770189432771018826227)

# ------------------------------------------------------------------
# Orders of the two Jacobians above K
# ------------------------------------------------------------------
qo = 1809251394333065553571917326471206521441306174399683558571672623546356726339
qt = 1809251394333065553414675955050290598923508843635941313077767297801179626051

K = SquaredKummerSurface(_zero, rosenhain=(_lam, _mu, _nu),
                          jacobian_order=2*qo, twist_order=2*qt)

# ------------------------------------------------------------------
# Generator point  (xDBL applied once to lift off special locus)
# ------------------------------------------------------------------
_gen_raw = (F(1), F(1), F(1), F(78525529738642755703105688163803666634))
gen = K.point(_gen_raw, validate=False).xDBL()

# ------------------------------------------------------------------
# Pre-generated random test points (as SquaredKummerPoints on K)
# ------------------------------------------------------------------
def _make_point(raw):
    return K.point(tuple(F(x) for x in raw))

random_points = [_make_point(raw) for raw in _raw_points]

# ------------------------------------------------------------------
# Convenience: canonical 2-torsion basis
# ------------------------------------------------------------------
basis = K.two_torsion_basis()


# ------------------------------------------------------------------
# Quick smoke test (run as script)
# ------------------------------------------------------------------
if __name__ == "__main__":
    from random import choice

    print("Null point on Kummer:", K.null_point().on_kummer())
    print("Generator on Kummer: ", gen.on_kummer())

    # Scalar multiplication sanity check: [qo] * gen == null point
    assert (qo * gen).is_zero(), "Order check failed"
    print("Order check passed: qo * gen == 0")

    # Profile of a doubled random point should be trivial (all True)
    P = choice(random_points)
    trivial = [True, True, True, True]
    prof = P.profile(basis)
    print(f"Profile of  P: {prof}  (trivial = {prof == trivial})")
    prof = P.xDBL().profile(basis)
    print(f"Profile of 2P: {prof}  (trivial = {prof == trivial})")
