## GOAL:
I'm building a class here for doing generic arithmetical things on (squared) Kummer surfaces.
If you are wondering why I'm not using "ordinary" theta structures, but squared theta structures:
 - this class includes Tate pairing computations using cubical arithmetic
 - for now, this requires working on squared Kummer surfaces

One nice application of this class is that it makes the computation of the Tate profile very easy,
which has ample applications when doing general arithmetic on Kummer surfaces, e.g. in sampling bases and so on

## TODO:
 - there is a bug in the computation of the Tate profile, which makes it non-linear for now (although with little effect on the profile)
 - due to some constraints, it is now focused on supersingular Jacobians of order (p+1)^2, although it should be easy to generalize to most other Jacobians, I think
 - there is almost no functionality for translating points on Jacobians to their Kummer points, eventhough this should be easy to do by applying the literature
