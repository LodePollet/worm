# PBC with linear length L == 3: resolved

## Summary

With periodic boundary conditions along an axis of linear length 3 (in a
lattice of dimension >= 2, e.g. `square` with `Lx=4, Ly=3`, or the
equivalent for `cubic`, `triangular`, `honeycomb`), the simulation would
intermittently abort during `test_conf()` with an error such as:

```
# TEST_CONF : error in chronology, two interactions at exactly the same
time on site <i> links : <j> <k>   times : <t> <t>
```

This branch previously added a guard rejecting `L == 3` under PBC
outright. That guard has since been removed: the underlying cause was
identified and the `test_conf()` diagnostic was fixed. `L == 3` under PBC
now works normally in `square`, `cubic`, `triangular`, and `honeycomb`.

## Root cause

For a periodic ring of length 3, any two neighbors of a site (reached by
stepping +1 and -1 along that axis) are themselves also direct neighbors
of each other -- the three sites form a triangle in which every pair is
bonded. This is a purely topological fact of an `L=3` ring, not specific
to any one lattice type.

If the worm dwells at a fixed imaginary time (`worm_at_stop != 0`) and
wanders around the 3-site triangle (which it can do because every pair of
the three sites is bonded), it can leave behind genuine, distinct kinks on
more than one bond of the triangle, all at the exact same time. From a
given site's point of view, its `+`-direction neighbor and its
`-`-direction neighbor then both report the *same* kink time -- not
through any floating-point coincidence or corruption, but because that
single kink genuinely sits on the bond connecting those two neighbor
sites to each other.

Directly inspecting the diagram at the moment of failure (site 6, in a
`Lx=4, Ly=3` run) confirmed this: the two "simultaneous" interactions
were not independent floating-point draws that coincidentally agreed, but
a *single* physical kink, recorded once on each of its two endpoint
sites, exactly as every kink always is:

```
Site 10: ... time : 3.319406723200699  link : 2   ...
Site  2: ... time : 3.319406723200699  link : 10  ...
```

Site 10's kink links to site 2, and site 2's kink links back to site 10,
at bit-identical times, because this is the same hop event on the 10-2
bond viewed from both ends. That bond only exists because `Ly=3` makes
site 10 and site 2 -- both neighbors of site 6, in opposite directions --
also neighbors of each other.

`test_conf()`'s consistency check (`src/worm.cpp`) assumed that a site's
distinct neighbor directions always refer to disjoint parts of the
diagram, so two different neighbor associations of the same site could
never legitimately point at the same imaginary time. That assumption is
false for `L == 3`, and the check flagged the legitimate shared-kink
case as corruption.

`dtol` and its scaling were tested up to 10000x with no effect on the
crash, confirming from the start that this was never a floating-point
tolerance problem.

### Why this was much rarer in 1D

The same topological degeneracy exists identically for a 1D `chain` of
length 3 (it is the same triangle graph, just with coordination number 2
instead of 4). Empirically, thousands of `test_conf()` checks per
micro-update ran clean on `chain, Lx=3` where the equivalent `square,
Ly=3` configuration failed within roughly 1000 updates. The likely reason
is dynamical: `MOVEWORM` picks its next target by comparing the
time-to-next-kink across *all* of a site's neighbor directions at once.
With `zcmax=4` (square), 4 directions compete every step, giving the worm
many more chances per step to wander off its current axis and
re-approach the same 3-site ring from a different direction. With
`zcmax=2` (chain), there are only 2 competing directions, so the
triangle-closing scenario is reached far less often, even though it was
never structurally impossible.

## The fix

`test_conf()` now distinguishes a genuine chronology error from a
legitimate shared-kink coincidence: when two neighbor associations
`assoc(j)` and `assoc(k)` of a site report the same time, it checks
whether the two neighbor sites `nb[i][j]` and `nb[i][k]` are themselves
bonded to each other and whether the two associations are literally the
mirror halves of that one bond's kink (`assoc(j)->link() == nb[i][k]` and
vice versa). Only if that is *not* the case is it flagged as corruption.
This is a fix to the diagnostic only -- `find_assoc_insert`,
`find_assoc_delete`, and the rest of the update logic were not modified.

## Validation

Before concluding the diagnostic fix was sufficient (as opposed to a
deeper bug in the update/sampling logic that the diagnostic had merely
been masking), the physics itself was validated independently:

- The exact configuration that previously crashed within ~1000 updates
  (`square`, `Lx=4, Ly=3`, hardcore bosons, mu=0, beta=4) was run for
  ~28M sweeps with a checked build (`test_conf()` running after every
  micro-update, throughout) with zero chronology failures.
- The resulting `G(k, i omega_n)` (48 points) and binned `G(k=0,tau)` (20
  bins) were compared against an independent exact-diagonalization
  calculation of the same system (full, untruncated Hilbert space, no
  windowing): 100% of points agreed within 3 sigma in both the frequency
  and imaginary-time domains.

This confirms the worm algorithm's sampling was correct at `L == 3` all
along -- the crash was purely a false positive in the consistency check,
not a sign of biased physics.
