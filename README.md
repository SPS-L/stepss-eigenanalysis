# STEPSS Eigenanalysis Reference Data

Reference spectra and the validation suite for small-signal stability analysis
in [STEPSS](https://stepss.sps-lab.org/).

**The analysis itself lives in RAMSES.** The engine reduces the linearised
differential-algebraic model to a state matrix, solves the eigenproblem, and
writes eigenvalues, damping ratios, participation factors and mode shapes,
triggered by an `EIG` disturbance or the `run_ssa` C entry. See the
[eigenanalysis user guide](https://stepss.sps-lab.org/user-guide/eigenanalysis/).

This repository holds the independently captured reference data that the engine
is checked against, plus the tests that check it. Nothing here is required to
run an analysis.

## Why a separate repository

The reference spectra were **not** produced by RAMSES. They were captured from an
independent implementation, so comparing the engine against them is a real
check rather than a self-fulfilling one. Keeping them outside the engine
repository also means these tests need neither a RAMSES licence nor the engine
itself: they run on numpy and pytest alone.

```sh
pip install numpy pytest
python -m pytest tests/ -v
```

## Layout

| Path | Contents |
|---|---|
| `fixtures/` | Exported Jacobians: four self-contained cases |
| `golden/` | Reference outputs for each case |
| `tests/` | The pytest suite |
| `capture_golden.py` | Regenerates `golden/` from `fixtures/` |

### Cases

| Case | States | Notes |
|---|---|---|
| `kundur_pss` | 70 | Kundur two-area with power system stabilisers |
| `kundur_nopss` | 70 | Identical but with `KSTAB = 0`, giving an unstable inter-area mode |
| `test` | 312 | Larger case, exercises scale |
| `1link_island` | 24 | Small case with a two-port and an island |

`example/py_*.dat` is **quarantined and deliberately excluded**: it holds 3,819
NaN entries in three equations of every exciter, and the capture script refuses
it.

### Reference outputs

For each case, `golden/` holds `<name>_Asys.txt` (the state matrix),
`<name>_eigs.txt` (eigenvalues), `<name>_pf.txt` (participation factors) and
`<name>_states.txt` (the state labels that index the participation rows).

Capturing `A_sys` separately from the eigenvalues is deliberate: it separates a
reduction error from an eigensolve error, which are otherwise indistinguishable
from a single failing comparison.

**Eigenvectors are not captured.** Their scale and phase are arbitrary, and for
a repeated eigenvalue they are not unique at all, so they cannot serve as a
reference.

## What the tests assert

`tests/test_golden.py` checks properties that hold on the reference data alone:
that the captured eigenvalues really are the spectrum of the captured state
matrix, that participation columns are normalised, and that the state labels
line up with the participation rows.

The assertion worth knowing about is the physical one. In the Kundur two-area
system the inter-area mode near 0.62 Hz is **unstable without the stabilisers**
and well damped with them:

| | inter-area | area 1 local | area 2 local |
|---|---|---|---|
| `kundur_nopss` | 0.625 Hz, zeta = **-0.0233** | 1.085 Hz, zeta = 0.099 | 1.116 Hz, zeta = 0.097 |
| `kundur_pss` | 0.624 Hz, zeta = **+0.1087** | 1.242 Hz, zeta = 0.288 | 1.295 Hz, zeta = 0.287 |

The sign flip is what a numerical regression cannot satisfy by accident, which
is why it is the assertion the engine's own release gate turns on. It reproduces
Kundur, *Power System Stability and Control*, Example 12.6.

## Degenerate eigenvalues

Identical device models with identical parameters produce identical poles, so
these spectra are heavily degenerate: 20 of the 70 modes in `kundur_nopss`, and
44 of the 312 in `test`.

**In a degenerate eigenspace the eigenvectors are not unique**, so per-mode
participation factors there are basis-dependent and differ legitimately between
implementations. Any comparison against `golden/*_pf.txt` must therefore be
restricted to modes whose gap to every other eigenvalue exceeds a tolerance.
Measured on `test`, the difference is 1.2e-08 across the 268 simple modes
against 7.6e-01 across the 44 degenerate ones, so a naive all-modes comparison
fails on correct code.

The engine reports this directly: `<name>_modes.dat` carries a simplicity flag
per mode.

## Regenerating the reference data

Only needed if a fixture changes:

```sh
for case in kundur_pss kundur_nopss test 1link_island; do
    python capture_golden.py fixtures/$case golden
done
```

`capture_golden.py` refuses to emit a reference whose eigenpair residual is too
large, so a silent pass is a real pass.

It shares no code with the Fortran engine it validates, which is what makes
agreement between them evidence rather than a tautology. It is a port of the
retired MATLAB implementation that originally produced these files, and
reproduces them to 1e-15 on the state matrix and 1e-14 on the eigenvalues. The
MATLAB sources themselves are gone; recover them from the history before commit
`1b604db` if you ever need to re-check that equivalence.

## License

Apache License 2.0. See [LICENSE](LICENSE). Copyright © Petros Aristidou.
[NOTICE](NOTICE) describes the licensing of proprietary components of the wider
STEPSS suite, which are not included here.

Developed and maintained by the
[Sustainable Power Systems Laboratory (SPS-L)](https://sps-lab.org/) at the
Cyprus University of Technology, under the direction of Dr. Petros Aristidou.

The Kundur fixtures derive from
[SPS-L/stepss-test-systems](https://github.com/SPS-L/stepss-test-systems)
(Apache-2.0). Please cite the original source of the system data:

> P. Kundur, *Power System Stability and Control*, McGraw-Hill, 1994
> (two-area system, Example 12.6).

For questions or support, please contact info@sps-lab.org.
