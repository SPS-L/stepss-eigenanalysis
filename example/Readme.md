# Archived example data

These files are kept as historical inputs and are **not** the way to run a
small-signal analysis.

RAMSES performs the analysis itself. Schedule an `EIG` disturbance and it writes
the modes, participation factors and mode shapes directly:

```python
ram.execSim(case, 0.0)
ram.addDisturb(0.001, "EIG 'ssa'")
ram.contSim(0.01)
```

See the [eigenanalysis user guide](https://stepss.sps-lab.org/user-guide/eigenanalysis/)
for the output format and the required settings, and the annotated notebook
under `examples/eigenanalysis/` in the `stepss` package for a full walkthrough
on the Kundur two-area system.

## About these files

`simply_load_and_run.ipynb` extracts a Jacobian into the four `*_val.dat`,
`*_eqs.dat`, `*_var.dat` and `*_struc.dat` files using the `JAC` disturbance.
That export still exists and is still useful if you want to drive your own
solver, but `getJac()` returns the same matrices directly as SciPy sparse
objects without the intermediate files.

The `test_*` and `1link_island_*` sets here are the originals of the copies
under `fixtures/`, which are what the test suite actually uses.

`py_*` is **defective and must not be used**: it holds 3,819 NaN entries in
three equations of every exciter in the case. `capture_golden.py` refuses it,
where the retired implementation would have consumed it silently and returned an
all-NaN spectrum.

## Notes that still apply

- The Jacobian export requires the synchronous reference frame
  (`$OMEGA_REF SYN ;`). Under the centre-of-inertia frame the export is skipped.
- The four-file export comes from the decomposed scheme (`$SCHEME DE`); the
  integrated scheme writes three, omitting `*_struc.dat`.
