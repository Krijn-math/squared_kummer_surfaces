"""
    Example: computing Tate profiles of degree 67 to check if R is in K(Fp) \ [67] K(Fp)
    all on the Kummer surface! We are lucky here, the existing sign ambiguity is no problem
    as we only need non-triviality
"""

from sage.all import GF
from kummer_surface import SquaredKummerSurface
from kummer_point import SquaredKummerPoint
from points import randompoints as _raw_points

# ------------------------------------------------------------------
# Base field
# ------------------------------------------------------------------
p = 2**246*3*67 - 1
ell = 67
F = GF(p)

# ------------------------------------------------------------------
# Kummer surface: null point and Rosenhain invariants
# ------------------------------------------------------------------
from_magma = [14563619561439671274064267294669858150005709431588288858052580823781116782630, 13454806426907545367845829201064487557983661848799930856637863004300863538119, 1, 1]
_zero = tuple([F(x) for x in from_magma])


# ------------------------------------------------------------------
# Orders of the two Jacobians above K up to doubling
# ------------------------------------------------------------------
qo = p+1
qt = p+1

K = SquaredKummerSurface(_zero, 
                        #  rosenhain=(_lam, _mu, _nu),
                          jacobian_order=qo, twist_order=qt)


if __name__ == "__main__":
        
    for i in range(10):
        R = K.random_point()
        if R.profile(ell).is_trivial():
            print( (((p+1) // ell) * R) == K.zero() )
        else:
            print( (((p+1) // ell) * R) != K.zero() )
            print( (ell*R).profile(ell).is_trivial())
        
        print("")
            
    