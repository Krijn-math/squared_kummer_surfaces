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
p = 2**246*3*67 - 1
ell = 67
F = GF(p)
R.<x> = F[]

# ------------------------------------------------------------------
# Kummer surface: null point and Rosenhain invariants
# ------------------------------------------------------------------
from_magma = [14563619561439671274064267294669858150005709431588288858052580823781116782630, 13454806426907545367845829201064487557983661848799930856637863004300863538119, 1, 1]
_zero = tuple([F(x) for x in from_magma])

qo = p+1
qt = p+1

# K = SquaredKummerSurface(_zero, jacobian_order=qo, twist_order=qt)

λ = F(12938376701500089469123999060444560237380085127671348976010005071974881035462)
μ = F(4618610882062834088149724968147708876591748255560510786118061150909387779558)
ν = F(1073370653850651346651397166382546615373880359440809096202890114345603658730)
ros = (λ, μ, ν)

K = SquaredKummerSurface(ros, jacobian_order=qo, twist_order=qt)


def rosenhain_to_hyperelliptic(ros):
    f = x*(x-1)*(x-ros[0])*(x-ros[1])*(x-ros[2])
    return HyperellipticCurve(f), f

def ros_to_jac(ros):
    H, f = rosenhain_to_hyperelliptic(ros)
    return H.jacobian(), f

# ------------------------------------------------------------------
# Quick smoke test (run as script)
# ------------------------------------------------------------------
if __name__ == "__main__":
    from random import choice

    H, f = rosenhain_to_hyperelliptic(ros)
    J = H.jacobian()
    J._ros = ros
            
    def random_curve_point(curve):
        while True:
            try:
                return H.lift_x(F.random_element())
                break
            except:
                pass
            
    def random_jacobian_element(jac):
        P1 = random_curve_point(jac.curve())
        P2 = random_curve_point(jac.curve())
        return J(P1) + J(P2)
    
    P = ((p+1) // ell)*random_jacobian_element(J)
    Q = random_jacobian_element(J)
    D = K(P - Q)
    P = K(P)
    Q = K(Q)