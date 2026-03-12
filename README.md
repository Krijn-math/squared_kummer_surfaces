## GOAL:
I'm building a class here for doing generic arithmetical things on (squared) Kummer surfaces.
If you are wondering why I'm not using "ordinary" theta structures, but squared theta structures:
 - this class includes Tate pairing computations using cubical arithmetic
 - for now, this requires working on squared Kummer surfaces

One nice application of this class is that it makes the computation of the Tate profile very easy,
which has ample applications when doing general arithmetic on Kummer surfaces, e.g. in sampling bases and so on.

One thing to be careful with when computing directly on the Kummer, e.g. not deriving from the Jacobian, is that we compute squares of Tate pairings, and can't differentiate between the Tate pairing or its inverse

## TODO:
 - when sampling points on the Kummer, we must (arbitrarily) choose a point difference in the Tate pairing computation, which effectively means we arbitrarily compute either t(P, R) or t(P, -R) = 1/t(P, R). When we derive these values from the Jacobian, we can resolve this issue
 - due to some constraints, it is now focused on supersingular Jacobians of order (p+1)^2, although it should be easy to generalize to most other Jacobians, I think
 - there is almost no functionality for translating points on Jacobians to their Kummer points, eventhough this should be easy to do by applying the literature
