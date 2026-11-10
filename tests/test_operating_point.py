"""Checks of the initialisation that need no MATLAB results."""

import numpy as np
import pytest

import wamces
from wamces.dynamics import Dynamics


@pytest.fixture(scope="module")
def undisturbed():
    return wamces.Model(disturbance_mw=0.0)


def test_sizes(undisturbed):
    m = undisturbed
    assert (m.nbus, m.nsm, m.ngfm, m.ngfl, m.ninj) == (3809, 618, 300, 500, 1418)
    assert (len(m.index_bs), len(m.index_pmu), len(m.index_pmu_bs)) == (2391, 1573, 1030)


def test_rest_without_disturbance(undisturbed):
    """Without the load step the model starts at rest."""
    dyn = Dynamics(undisturbed)
    dx, mem_in, sig = dyn.evaluate(dyn.pack(undisturbed.x0), undisturbed.m0)
    assert np.max(np.abs(dx)) < 1e-10
    assert abs(sig["dfsys"]) < 1e-15 and abs(sig["dfw"]) < 1e-15
    for k in mem_in:
        assert np.max(np.abs(mem_in[k] - undisturbed.m0[k])) < 1e-10


def test_voltages_match_power_flow(undisturbed):
    """At the start the network gives the voltages of the power flow."""
    m = undisturbed
    dyn = Dynamics(m)
    _, _, sig = dyn.evaluate(dyn.pack(m.x0), m.m0)
    buses = np.concatenate([m.index_inj, m.index_pmu_bs])
    V_pf = m.Vm[buses] * np.exp(1j * m.Va[buses])
    assert np.max(np.abs(sig["Vbus"][buses] - V_pf)) < 1e-10


def test_disturbance_lowers_frequency():
    r = wamces.Model().simulate(t_end=2.0, save_every=100)
    assert r.f_coi[-1] < wamces.model.FN - 1e-3


def test_initialise_twice_gives_the_same_case(undisturbed):
    """A second initialise starts from the ratings of the data again."""
    m = undisturbed
    sr, x0, n_over = m.sm.Sr.copy(), m.x0["sm_eq1"].copy(), len(m.overloaded)
    m.initialise()
    assert n_over == 39 and len(m.overloaded) == n_over
    assert np.array_equal(m.sm.Sr, sr) and np.array_equal(m.x0["sm_eq1"], x0)
