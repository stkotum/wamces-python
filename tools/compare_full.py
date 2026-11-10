"""Compare every state at every step with a full Simulink log.

Usage:
    python tools/compare_full.py states_published.mat published
    python tools/compare_full.py states_variant.mat variant

The files ``states_<case>.mat`` come from ``tools/make_reference.m`` with WAMCES_FULL=1.
The script needs h5py. It prints one Markdown table row per state group. Each row gives
the largest absolute error and that error relative to the largest value of the group.
The published case needs about 2 GB of memory.
"""

import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))

import wamces  # noqa: E402
from wamces.dynamics import SIMULINK_BLOCKS  # noqa: E402
from test_simulink import variant_model  # noqa: E402


def main(path, case):
    f = h5py.File(path, "r")
    blocks = ["".join(chr(c) for c in f[r][()].ravel()) for r in f["block"][()].ravel()]
    t_ref = f["t"][()].ravel()
    model = wamces.Model() if case == "published" else variant_model()
    result = model.simulate(t_end=float(t_ref[-1]), save_every=1)
    assert np.allclose(result.t, t_ref, rtol=0, atol=1e-9)
    print("| Group | States | Largest error | Relative error |")
    print("|---|---|---|---|")
    worst = 0.0
    for k, block in enumerate(blocks, start=1):
        name = SIMULINK_BLOCKS[block.replace("wamces_model/", "", 1)]
        ref = f["s%03d" % k][()].T  # MATLAB v7.3 stores the transpose
        err = np.max(np.abs(result.x[name] - ref))
        rel = err / max(np.max(np.abs(ref)), 1e-300)
        worst = max(worst, rel)
        print("| `%s` | %d | %.1e | %.1e |" % (name, ref.shape[1], err, rel))
    print("\nworst relative error %.1e over %d steps" % (worst, len(t_ref)))


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[2] not in ("published", "variant"):
        sys.exit(__doc__)
    main(*sys.argv[1:])
