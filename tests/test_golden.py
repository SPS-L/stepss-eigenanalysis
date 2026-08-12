"""Checks that hold on the goldens alone, with no RAMSES and no MATLAB."""
import numpy as np, pathlib, pytest

GOLDEN = pathlib.Path(__file__).parent.parent / "golden"
CASES = ["kundur_pss", "kundur_nopss", "test", "1link_island"]


def load_asys(name):
    t = np.loadtxt(GOLDEN / f"{name}_Asys.txt")
    n = int(max(t[:, 0].max(), t[:, 1].max()))
    a = np.zeros((n, n))
    a[t[:, 0].astype(int) - 1, t[:, 1].astype(int) - 1] = t[:, 2]
    return a


def load_eigs(name):
    e = np.loadtxt(GOLDEN / f"{name}_eigs.txt")
    return e[:, 0] + 1j * e[:, 1]


@pytest.mark.parametrize("name", CASES)
def test_eigenvalues_solve_the_golden_state_matrix(name):
    """The captured spectrum must be the spectrum of the captured matrix."""
    a, w_gold = load_asys(name), load_eigs(name)
    w = np.linalg.eigvals(a)
    order = lambda z: np.lexsort((z.imag, z.real))
    err = np.abs(w[order(w)] - w_gold[order(w_gold)]).max()
    assert err / np.abs(w_gold).max() < 1e-12, f"{name}: {err:.3e}"


@pytest.mark.parametrize("name", CASES)
def test_participation_columns_peak_at_one(name):
    p = np.loadtxt(GOLDEN / f"{name}_pf.txt")
    assert np.allclose(p.max(axis=0), 1.0)
    assert (p >= 0).all()


@pytest.mark.parametrize("name", CASES)
def test_state_labels_match_participation_rows(name):
    p = np.loadtxt(GOLDEN / f"{name}_pf.txt")
    labels = (GOLDEN / f"{name}_states.txt").read_text().splitlines()
    assert len(labels) == p.shape[0]


def test_pss_flips_the_interarea_damping_sign():
    """Kundur Example 12.6: the inter-area mode is unstable without the PSS.

    This is the assertion that a numerical regression cannot satisfy by
    accident, so it is the one worth stating loudest.
    """
    def interarea(name):
        w = load_eigs(name)
        f = np.abs(w.imag) / (2 * np.pi)
        band = w[(f > 0.4) & (f < 0.9) & (w.imag > 0)]
        assert len(band) == 1, f"{name}: expected one inter-area mode, got {len(band)}"
        return -band[0].real / abs(band[0])

    assert interarea("kundur_nopss") < 0, "no-PSS inter-area mode should be unstable"
    assert interarea("kundur_pss") > 0.05, "PSS should damp the inter-area mode"
