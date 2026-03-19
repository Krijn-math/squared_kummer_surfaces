## GOAL:
I'm building a class here for doing generic arithmetical things on Jacobians and (squared) Kummer surfaces.
On both objects, we are using cubical arithmetic to compute pairings, e.g. Tate, Weil, perhaps cokernel pairings.

If you are wondering why I'm not using "ordinary" theta structures, but squared theta structures:
 - this class includes Tate pairing computations using cubical arithmetic
 - for now, this requires working on squared Kummer surfaces

One nice application of this class is that it makes the computation of the Tate profile very easy,
which has ample applications when doing general arithmetic on Kummer surfaces, e.g. in sampling bases and so on.

One thing to be careful with when computing directly on the Kummer, e.g. not deriving from the Jacobian, is that we compute squares of Tate pairings, and can't differentiate between the Tate pairing or its inverse

## EXAMPLES:
There are a few examples on computing pairings on jacobians and kummers, and one for computing profiles. 
See the files that start with `example_` for more explanation.
Furthermore, for my note on the Tate pairing (ePrint 2025/477), there is a folder with examples for most use cases.

## TODO:
 - due to some constraints, it is now focused on supersingular Jacobians of order (p+1)^2, although it should be easy to generalize to most other Jacobian
 - there is only limited functionality for translating points on Jacobians to their Kummer points; we rely on rosenhain invariants to make this work
 - - hence, usually, you should precompose your Jacobian with an isomorphism to one in Rosenhain form
