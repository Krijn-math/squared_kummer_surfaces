# TODO

## Linearity of Tate pairing
Current approach doesnt seem to respect the linearity, not necessary for profile but would be nice to have

## Quick check on multiples
When generating a basis, it is now difficult to check that P2 is not a multiple of P1.
Can we somehow improve this with Tate pairings and dlogs?

# Backlog

## 1. Derive Rosenhain invariants from null point
Implemented as `K.rosenhain_invariants(bit=0)`. Given null point (a,b,c,d) and dual (A,B,C,D) = hadamard(a,b,c,d): compute alpha = sqrt(C*D/(A*B)), set ratio = (1+alpha)/(1-alpha), then lam = ac/(bd), mu = c*ratio/d, nu = a*ratio/b. One square root; bit=1 negates alpha to select the conjugate triple.

## 2. Derive Mij addition matrices from null point alone
Currently requires Rosenhain invariants (computed via formulas in `kummer_surface.py` `_compute_addition_matrices`).
Once backlog item 1 is done, Mij can be computed from zero alone.

## 3. Verify M[3][5] delta entry
In `kummer_surface.py` `_compute_addition_matrices`, the `delta` variable for `M[3][5]` is:

    delta = r(3,6)/r(3,5) * m[2]/m[4]

i.e. `(w[3]-w[6])/(w[3]-w[5]) * m[2]/m[4]`. This matches `trying_out_matrices.py` row 2 col 0, but should be verified against the original Magma derivation since `mus[1]` appears in some entries of `M[3][5]` which has an unusual structure compared to the other matrices.

## 4. RotK pairing method
The `helper.py` / `origin.py` approach using Mumford coordinates `(u0, u1, v0^2)` via Rosenhain invariants. Currently removed; may be added later as `P.profile(basis, method='rotk')`.

## 5. Optimized cubical pairing
`optcubical.py` had a more efficient profile using basis `(1,2),(2,4),(2,3),(4,5)` with precomputed lambdas and simplified Legendre checks. Can be re-added as an optimization layer.

## 6. Operation counting / benchmarking
The old `finite_field.py` tracked mul/sqr/inv/leg counts. Should be re-added as an optional instrumentation layer on top of Sage field ops (e.g. via a wrapper field class or context manager).
