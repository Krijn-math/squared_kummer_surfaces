# TODO


## Quick check on multiples
When generating a basis, it is now difficult to check that P2 is not a multiple of P1.
Can we somehow improve this with Tate pairings and dlogs?

# Backlog

## 3. Verify M[3][5] delta entry
In `kummer_surface.py` `_compute_addition_matrices`, the `delta` variable for `M[3][5]` is:

    delta = r(3,6)/r(3,5) * m[2]/m[4]

i.e. `(w[3]-w[6])/(w[3]-w[5]) * m[2]/m[4]`. This matches `trying_out_matrices.py` row 2 col 0, but should be verified against the original Magma derivation since `mus[1]` appears in some entries of `M[3][5]` which has an unusual structure compared to the other matrices.

## 5. Optimized cubical pairing
`optcubical.py` had a more efficient profile using basis `(1,2),(2,4),(2,3),(4,5)` with precomputed lambdas and simplified Legendre checks. Can be re-added as an optimization layer.

## 6. Operation counting / benchmarking
The old `finite_field.py` tracked mul/sqr/inv/leg counts. Should be re-added as an optional instrumentation layer on top of Sage field ops (e.g. via a wrapper field class or context manager).
