#!/usr/bin/env python3
"""Capture reference outputs of the small-signal reduction as plain text.

Reads an exported Jacobian (``<prefix>_val.dat``, ``_eqs.dat``, ``_var.dat``),
reduces it to the state matrix, solves the dense eigenproblem, and writes the
reference files that ``tests/test_golden.py`` and the RAMSES release gates are
checked against.

This is deliberately an **independent** implementation. It shares no code with
the Fortran engine it validates, so agreement between the two is evidence
rather than a tautology. It is a direct port of the retired MATLAB reference
implementation and reproduces its output to better than 1e-14.

Usage::

    python capture_golden.py fixtures/kundur_pss golden
    python capture_golden.py fixtures/test golden

Outputs, all 1-based and full double precision:

``<name>_Asys.txt``
    ``i j value`` per nonzero of the state matrix. Captured separately from the
    eigenvalues on purpose: it separates a reduction error from an eigensolve
    error, which one comparison alone cannot distinguish.
``<name>_eigs.txt``
    ``real imag`` per eigenvalue.
``<name>_pf.txt``
    Dense participation factors, rows = differential states, columns = modes,
    each column normalised to peak at 1.
``<name>_states.txt``
    ``index family device variable`` per differential state, giving the row
    labels of ``_pf.txt``.

Modes are sorted by descending real part, then descending imaginary part, so the
files are reproducible.

Eigenvectors are **not** captured. Their scale and phase are arbitrary, and for
a repeated eigenvalue they are not unique at all, so they cannot serve as a
reference. Participation factors are only comparable for simple eigenvalues for
the same reason: identical device models produce identical poles, and these
systems are heavily degenerate.
"""

import os
import sys

import numpy as np
from scipy.linalg import eig
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


def build_state_matrix(prefix):
    """Return the state matrix, the differential state indices, and their labels.

    The linearised model is ``dx/dt = fx dx + fy dy`` and ``0 = gx dx + gy dy``.
    Eliminating the algebraic variables gives the Schur complement
    ``A_sys = fx - fy gy^-1 gx``, which exists only if ``gy`` is nonsingular.
    """
    val = np.loadtxt(prefix + "_val.dat")
    if val.ndim != 2 or val.shape[1] < 3:
        raise SystemExit(
            "%s_val.dat is not a plain 3-column coordinate file. Non-numeric "
            "entries such as NaN in the value column are the usual cause." % prefix)
    if not np.isfinite(val[:, 2]).all():
        raise SystemExit(
            "%s_val.dat holds %d non-finite entries (NaN or Inf) of %d; the "
            "spectrum would be all NaN."
            % (prefix, int((~np.isfinite(val[:, 2])).sum()), val.shape[0]))

    n = int(max(val[:, 0].max(), val[:, 1].max()))
    # coo -> csr sums duplicate entries, which the assembly does produce.
    S = coo_matrix((val[:, 2], (val[:, 0].astype(int) - 1, val[:, 1].astype(int) - 1)),
                   shape=(n, n)).tocsr()

    gamma = np.loadtxt(prefix + "_eqs.dat", usecols=5, dtype=int)
    names = np.loadtxt(prefix + "_var.dat", usecols=(1, 2, 3, 4), dtype=str)

    # Permute rows so equation i differentiates state i, making E diagonal. The
    # algebraic rows land in whatever order remains, which does not matter: gy
    # and gx receive the same permutation and (P gy)^-1 (P gx) = gy^-1 gx.
    perm = np.empty(n, dtype=int)
    claimed = np.zeros(n, dtype=bool)
    for j, g in enumerate(gamma):
        if g != 0:
            perm[g - 1] = j
            claimed[g - 1] = True
    free = [j for j, g in enumerate(gamma) if g == 0]
    for slot, j in zip(np.flatnonzero(~claimed), free):
        perm[slot] = j
    S = S[perm, :]

    dif = np.flatnonzero(names[:, 0] == "d")
    alg = np.flatnonzero(names[:, 0] != "d")

    fx, fy = S[dif][:, dif], S[dif][:, alg]
    gx, gy = S[alg][:, dif], S[alg][:, alg]

    # Factor gy once and back-substitute, rather than re-solving per column.
    gygx = splu(gy.tocsc()).solve(gx.toarray())
    a_sys = np.asarray((fx - fy @ gygx).todense()) if hasattr(fx - fy @ gygx, "todense") \
        else np.asarray(fx - fy @ gygx)
    return a_sys, dif, names


def capture(prefix, outdir):
    name = os.path.basename(prefix)
    a_sys, dif, names = build_state_matrix(prefix)
    nx = a_sys.shape[0]

    lam, vl, vr = eig(a_sys, left=True, right=True)

    order = np.lexsort((-lam.imag, -lam.real))
    lam, vl, vr = lam[order], vl[:, order], vr[:, order]

    # Participation, normalised per mode so the arbitrary eigenvector scaling
    # cancels. abs() is taken before normalising, matching the reference.
    pf = np.abs(np.conj(vl) * vr)
    pf /= pf.max(axis=0)

    # Independent of any reference, so it must hold here.
    res = np.linalg.norm(a_sys @ vr - vr * lam, "fro") / np.linalg.norm(a_sys, "fro")
    print("%s: Nx=%d, relative eigenpair residual %.3e" % (name, nx, res))
    if res > 1e-10:
        raise SystemExit(
            "Eigenpair residual %.3e is too large to use as a reference." % res)

    os.makedirs(outdir, exist_ok=True)
    j, i = np.meshgrid(np.arange(nx), np.arange(nx))
    nz = a_sys != 0.0
    np.savetxt(os.path.join(outdir, name + "_Asys.txt"),
               np.column_stack([i[nz] + 1, j[nz] + 1, a_sys[nz]]),
               fmt=["%d", "%d", "%.17e"])
    np.savetxt(os.path.join(outdir, name + "_eigs.txt"),
               np.column_stack([lam.real, lam.imag]), fmt="%.17e")
    np.savetxt(os.path.join(outdir, name + "_pf.txt"), pf, fmt="%.17e")
    with open(os.path.join(outdir, name + "_states.txt"), "w") as fh:
        for k in dif:
            fh.write("%d %s %s %s\n" % (k + 1, names[k, 1], names[k, 2], names[k, 3]))

    print("%s: wrote 4 golden files to %s" % (name, outdir))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    capture(sys.argv[1], sys.argv[2])
