# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Building

All build artifacts go in `build/`. The source is in `src/`.

```bash
cd build
cmake ../src -DLATTICE=chain        # required: chain | square | cubic | ladder | triangular | honeycomb
make -j4
```

Key CMake flags (all optional):

| Flag | Default | Effect |
|---|---|---|
| `-DLATTICE=<name>` | `chain` | Lattice geometry — **must match the physics you want** |
| `-DDSTRUC=AVL\|LIST\|LIST_STACK` | `LIST` | Internal data structure for the diagram |
| `-DDEBUG=ON` | `OFF` | Enables `DEBUGMODE`: runs `test_conf()` after every update |
| `-DUNIFORM=OFF` | `ON` | Uniform (`UNISYS`) vs site/bond-inhomogeneous parameters |
| `-DCWINDOW=ON` | `OFF` | Canonical window for Green's function measurement |

This produces two executables: `qmc_worm` (single-core) and `qmc_worm_mpi` (MPI).

## Running

```bash
./build/qmc_worm parameter_files/BoseHubbard.ini
mpiexec -n 4 ./build/qmc_worm_mpi parameter_files/BoseHubbard.ini
```

Resume from checkpoint: `./qmc_worm job.clone.h5`

Output is HDF5 (`job.out.h5`). Helper scripts in `tools/` extract results:

```bash
python tools/print_params.py job.out.h5
python tools/print_corrfun.py job.out.h5
python tools/print_update_statistics.py job.out.h5
```

## Tests

Benchmarks against exact Lanczos diagonalization live in `test_mpi/`:

```bash
cd test_mpi/BoseHubbard_GrandCanonical
bash test_serial.sh chain        # runs cmake + make + simulation, then compare.py checks results
```

Subdirectories: `BoseHubbard_GrandCanonical`, `BoseHubbard_Canonical`, `XXZ_SpinOneHalf`, `XXZ_SpinOne`.

## Architecture

### Entry points
- `src/worm.run.cpp` — single-core `main()`
- `src/worm.run_mpi.cpp` — MPI `main()`

Both construct a `worm` object and hand off to the ALPSCore MC driver loop.

### Core class: `worm` (`src/worm.hpp`, `src/worm.cpp`)
Inherits `alps::mcbase`. Owns the diagram (a `Diagram_type` container of `Element`s per site), the worm head/tail iterators, and all measurement accumulators. Key members:
- `worm_diag` — bool distinguishing the diagonal (no worm) from off-diagonal (worm present) sector
- `update_prob_cuml[]` — cumulative probabilities for selecting each update

### MC updates (`src/worm.update.cpp`)
Five updates: `INSERTWORM`, `GLUEWORM`, `MOVEWORM`, `INSERTKINK`, `DELETEKINK`. Called from `worm::update()` which dispatches based on `worm_diag`.

### Data structure (`src/worm.Element.hpp`)
`Element` stores a single vertex: imaginary time, site index (`mLink`), occupation before/after (`mBefore`/`mAfter`), color (boson kink = +1, worm = −1, dummy = 0), and `mAssoc[ZC]` — iterators to the nearest neighbor elements just ahead in time on each neighboring site. Three interchangeable containers selectable at compile time via `-DDSTRUC`: `AVL_diagram`, `list_diagram`, `list_stack_diagram`.

The `operator<` tie-breaking for equal-time elements (comparing `after()` vs `before()`) is an intentional physics-based workaround: two worm elements at the same time on the same site form a chain via continuous occupancy, making the ordering unambiguous. Do not treat this as a bug.

### Models (`src/model.hpp`)
Abstract base `model` with two concrete implementations: `BoseHubbard` and `XXZ`. Each provides diagonal/off-diagonal site and bond weights. Compiled with `UNISYS` (uniform parameters, scalars) or without (inhomogeneous, vectors indexed by site/bond).

### Lattice (`src/lattice.hpp`)
Template `lattice<DIM, N_BASIS>` selected at compile time via `-DLATTICE=`. Provides neighbor maps, bond lists, and boundary condition logic. `LATTICE::zcmax` sets the compile-time maximum coordination number used to size `mAssoc[ZC]` in every `Element`.

### `src_xml/`
Older variant of the source — kept for reference but not the active build target.
