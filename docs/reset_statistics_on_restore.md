# Resetting statistics on restore: `--reset-statistics` CLI flag

## Symptom

The original mechanism for forcing a statistics reset when resuming from
a checkpoint was an `alps::params`-defined `reset_statistics` parameter,
checked inside `worm::load()`:

```cpp
reset_statistics = parameters["reset_statistics"];
if (reset_statistics == 1) {
  force_reset_statistics();
}
```

This could never actually be triggered. `worm::load()` is only called
when `parameters.is_restored()` is true (i.e. the simulation was
constructed by loading an HDF5 checkpoint archive as its first
argument), and `worm::define_parameters()` returns immediately in
exactly that case:

```cpp
void worm::define_parameters(parameters_type & parameters) {
    if (parameters.is_restored()) {
        return;
    }
    ...
}
```

Since `define_parameters()` never runs on restore, no parameter --
`reset_statistics` included -- can be freshly defined or overridden via
the command line or an INI file at that point. `parameters["reset_statistics"]`
reads back only whatever value happened to be stored in the checkpoint
archive from the *original* run that created it, which is always `0`
(its default), since nobody sets `reset_statistics=1` on the run that
produces the checkpoint in the first place. This was verified directly:

```
$ ./qmc_worm job.clone.h5 --reset_statistics=1
...
# reset statistics when restored (requires hack)    : 0
```

The CLI override is silently ignored regardless of the exact syntax
tried (`--reset_statistics=1`, `reset_statistics=1`, before or after the
checkpoint path). This is a general property of how `alps::params`
handles restored archives, not specific to this parameter -- so
`force_reset_statistics()` was effectively dead code, unreachable via
the documented resume workflow.

## The fix

Statistics reset-on-restore is now handled by a `--reset-statistics`
command-line flag, parsed manually in `main()` (`worm.run.cpp` and
`worm.run_mpi.cpp`) independently of `alps::params` -- this works
regardless of whether the run is a fresh start or a restore, since it
never goes through the parameter-definition machinery at all:

```cpp
bool reset_statistics_flag = std::any_of(argv + 1, argv + argc,
    [](const char* a) { return std::string(a) == "--reset-statistics"; });
...
if (parameters.is_restored()) {
    sim.load(checkpoint_file);
    if (reset_statistics_flag) {
        sim.force_reset_statistics();
    }
}
```

Usage: `./qmc_worm job.clone.h5 --reset-statistics`.

The `reset_statistics` member variable remains, but is now purely
informational (set by `force_reset_statistics()` itself, printed by
`print_params()`) rather than being read from `alps::params`.

## Validation

Confirmed by comparing measurement counts across three runs sharing one
checkpoint:

- Initial run: `Kinetic_Energy` count = 457243.
- Resume *without* `--reset-statistics`: count = 500000 (continues from
  457243, as expected -- no reset).
- Resume *with* `--reset-statistics`: prints "resetting statistics",
  `print_params()` shows `statistics were force-reset on restore : 1`,
  and the count (489353) reflects a genuinely fresh accumulation cycle
  rather than a continuation from 457243 -- consistent with
  `force_reset_statistics()` rewinding both the ALPS accumulators and
  the internal `sweeps` counter back to `thermalization_sweeps`.

Verified for both `qmc_worm` (single-core) and `qmc_worm_mpi`.
