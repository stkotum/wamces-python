"""Comparison with the published Simulink model.

The file tests/reference/wamces_reference.mat comes from tools/make_reference.m. It holds
all states at a few steps and three columns of every state block over time. The
Python states must agree with Simulink to 1e-10 of the largest value of each block.
The agreement found is 8e-12 or better (docs/validation.md).
"""

from pathlib import Path

import numpy as np
import pytest

import wamces
from wamces.dynamics import SIMULINK_BLOCKS

REFERENCE = Path(__file__).parent / "reference" / "wamces_reference.mat"
RTOL = 1e-10


@pytest.fixture(scope="module")
def reference():
    sio = pytest.importorskip("scipy.io")
    if not REFERENCE.exists():
        pytest.skip("no reference file, run tools/make_reference.m")
    return sio.loadmat(REFERENCE, squeeze_me=True, struct_as_record=False)


def variant_model():
    """The variant case of tools/make_reference.m."""
    m = wamces.Model(disturbance_bus=1, disturbance_mw=1000.0, kw=20.0)
    m.gfl.P0 = m.gfl.P0 + 0.1
    return m


def compare(R, result):
    """Largest error of each state group relative to the largest value of its block."""
    steps = np.atleast_1d(R.snap_steps).astype(int) // int(R.every)
    errors = {}
    for k, block in enumerate(np.atleast_1d(R.block), start=1):
        name = SIMULINK_BLOCKS[block.replace("wamces_model/", "", 1)]
        snap = np.atleast_2d(getattr(R, "snap%02d" % k))
        cols = np.atleast_1d(getattr(R, "cols%02d" % k)).astype(int) - 1
        trace = getattr(R, "trace%02d" % k).reshape(len(result.t), len(cols))
        py = result.x[name]
        scale = max(np.max(np.abs(snap)), np.max(np.abs(trace)), 1e-300)
        err = max(np.max(np.abs(py[steps] - snap)), np.max(np.abs(py[:, cols] - trace)))
        errors[name] = err / scale
    return errors


def check(errors):
    worst = max(errors, key=errors.get)
    assert errors[worst] < RTOL, "%s differs by %.3e" % (worst, errors[worst])


def test_variant(reference):
    R = reference["variant"]
    result = variant_model().simulate(t_end=10.0, save_every=int(R.every))
    assert np.allclose(result.t, R.t, rtol=0, atol=1e-9)
    check(compare(R, result))


@pytest.mark.slow
def test_published(reference):
    R = reference["published"]
    result = wamces.Model().simulate(t_end=50.0, save_every=int(R.every))
    assert np.allclose(result.t, R.t, rtol=0, atol=1e-9)
    check(compare(R, result))
