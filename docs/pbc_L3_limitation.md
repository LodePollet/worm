# PBC with linear length L == 3 is not supported

## Symptom

With periodic boundary conditions along an axis of linear length 3 (in a
lattice of dimension >= 2, e.g. `square` with `Lx=4, Ly=3`, or the
equivalent for `cubic`, `triangular`, `honeycomb`), the simulation would
intermittently abort during `test_conf()` with an error such as:

```
# TEST_CONF : error in chronology, two interactions at exactly the same
time on site <i> links : <j> <k>   times : <t> <t>
```

The lattice constructors now reject `L == 3` under PBC outright (in
addition to the pre-existing rejection of `L < 3`, which causes double
counting of bonds), rather than letting the simulation run and
occasionally hit this error.

## Root cause

For a periodic ring of length 3, any two neighbors of a site (reached by
stepping +1 and -1 along that axis) are themselves also direct neighbors
of each other -- the three sites form a triangle in which every pair is
bonded. This is a purely topological fact of an `L=3` ring and is *not*
specific to any one lattice type: it applies identically to a 1D `chain`
of length 3 and to any linear dimension of length 3 embedded in a
higher-dimensional lattice.

The diagram bookkeeping in `find_assoc_insert()` (`src/worm.cpp`) walks,
for a newly inserted kink on `cursite`, each neighbor direction `j`
independently to find the nearest existing kink in imaginary time on that
neighbor's own operator string, and to propagate the corresponding
"opposite direction" association back. This logic -- and the `test_conf()`
consistency check that later verifies it -- implicitly assumes that a
site's distinct neighbor directions refer to *disjoint* parts of the
diagram, so that two different neighbor associations of the same site can
never legitimately point at the same imaginary time.

That assumption is false whenever `L == 3`: if the worm dwells at a fixed
imaginary time (`worm_at_stop != 0`) and wanders around the 3-site
triangle (which it can do because every pair of the three sites is
bonded), it can leave behind genuine, distinct kinks on more than one
bond of the triangle, all at the exact same time. From `cursite`'s point
of view, its `+`-direction neighbor and its `-`-direction neighbor then
both report the *same* kink time -- not through any floating-point
coincidence or corruption, but because that single kink genuinely sits on
the bond connecting those two neighbor sites to each other (see below).

### Verifying it is not roundoff

Directly inspecting the diagram at the moment of failure (site 6, in a
`Lx=4, Ly=3` run) confirmed that the two "simultaneous" interactions are
not independent floating-point draws that coincidentally agree, but a
*single* physical kink, recorded once on each of its two endpoint sites,
exactly as every kink always is:

```
Site 10: ... time : 3.319406723200699  link : 2   ...
Site  2: ... time : 3.319406723200699  link : 10  ...
```

Site 10's kink links to site 2, and site 2's kink links back to site 10,
at bit-identical times, because this is the same hop event on the 10-2
bond viewed from both ends. That bond only exists because `Ly=3` makes
site 10 and site 2 -- both neighbors of site 6, in opposite directions --
also neighbors of each other. `find_assoc_insert` and `test_conf()` were
not written to expect that.

### Why this is much rarer in 1D

The same topological degeneracy exists identically for a 1D `chain` of
length 3 (it is the same triangle graph, just with coordination number 2
instead of 4). Empirically, though, thousands of `TEST_CONF` checks per
micro-update run clean on `chain, Lx=3` where the equivalent `square,
Ly=3` configuration fails within roughly 1000 updates.

The likely reason is dynamical, not structural: `MOVEWORM` picks its next
target by comparing the time-to-next-kink across *all* of a site's
neighbor directions at once (`src/worm.update.cpp`, the loop over
`zc[isite]` that selects `nb_next`). With `zcmax=4` (square), 4 directions
compete every step, giving the worm many more chances per step to wander
off its current axis and re-approach the same 3-site ring from a
different direction. With `zcmax=2` (chain), there are only 2 competing
directions -- forward and backward along the same line -- so the worm's
exploration is far more constrained and the triangle-closing scenario is
reached far less often, even though it is not structurally impossible.

## What would be required to support L == 3

Making `L == 3` correct (rather than merely rejecting it) would require
`find_assoc_insert`, `find_assoc_delete`, and the association-based
consistency check in `test_conf()` to stop assuming a site's neighbor
directions are diagram-disjoint, and instead handle the case where two
different neighbor associations legitimately alias the same kink. This is
a structural change to the association bookkeeping, not a matter of
adjusting a tolerance parameter (`dtol` and its scaling were tested up to
10000x with no effect, confirming the issue is not a floating-point
tolerance problem). Until that is done, `L == 3` under PBC is rejected at
lattice construction time in `square`, `cubic`, `triangular`, and
`honeycomb`.
