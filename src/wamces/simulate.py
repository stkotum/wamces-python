"""Time simulation with the fixed-step solver ode3 of Simulink.

Simulink picks ode3 for the published setting "FixedStepAuto" with a step of 0.01 s.
The method of Bogacki and Shampine takes three evaluations of the model per step.
The code follows the ode3 of Simulink, including the order of its sums.

Each step starts at a major time step. The model is evaluated there. This evaluation gives
the first stage. The two later stages are minor time steps.

A Memory block keeps the same output over the whole step. It gives the new value from
the next major time step on. This order reproduces the Simulink states to about 1e-12.
A Memory block that changes its output already in the minor time steps gives errors of
a few percent.
"""

import numpy as np

from .dynamics import Dynamics, MEMORY_GROUPS, STATE_GROUPS
from .model import FN

# Butcher tableau of ode3 as in Simulink.
_B = ((1.0 / 2.0, 0.0, 0.0), (0.0, 3.0 / 4.0, 0.0), (2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0))


class Result:
    """The result of a simulation.

    Attributes
    ----------
    t : array
        The times of the saved states in s.
    x : dict of arrays
        One array per saved state group, with one row per saved time.
    t_signals : array
        The times of the signals in s. The signals are saved at every step.
    dfsys : array
        Mean speed deviation of the synchronous machines in per unit.
    dfw : array
        Frequency deviation that the central unit measures, in per unit.
    """

    def __init__(self, t, x, t_signals, signals):
        self.t = t
        self.x = x
        self.t_signals = t_signals
        self.dfsys = signals["dfsys"]
        self.dfw = signals["dfw"]

    @property
    def f_sm(self):
        """Frequency of each synchronous machine in Hz."""
        return FN * (1.0 + self.x["sm_domega"])

    @property
    def f_pmu(self):
        """Frequency that each PMU measures in Hz."""
        return self.x["pmu_f"]

    @property
    def f_coi(self):
        """Mean frequency of the synchronous machines in Hz, at every step."""
        return FN * (1.0 + self.dfsys)

    @property
    def f_wide_area(self):
        """Frequency that the central unit measures in Hz, at every step."""
        return FN * (1.0 + self.dfw)


def simulate(model, t_end=50.0, h=0.01, save_every=1, save=None):
    """Simulate ``model`` from 0 to ``t_end`` with the fixed step ``h``.

    Parameters
    ----------
    model : wamces.Model
    t_end, h : float
        End time and step in s. The published model uses 50 s and 0.01 s.
    save_every : int
        Save the states at every ``save_every``-th step. The signals are saved at
        every step.
    save : list of str, optional
        The state groups to save. The default saves all groups of ``STATE_GROUPS``.
        All groups at every step of 50 s take about 840 MB.

    Returns
    -------
    Result
    """
    dyn = Dynamics(model)
    groups = list(STATE_GROUPS) if save is None else list(save)
    x = dyn.pack(model.x0)
    mem = {k: np.array(model.m0[k], copy=True) for k in MEMORY_GROUPS}
    nsteps = int(round(t_end / h))
    nsave = nsteps // save_every + 1
    out = {g: np.empty((nsave, dyn.layout[g].stop - dyn.layout[g].start)) for g in groups}
    sig = {"dfsys": np.empty(nsteps + 1), "dfw": np.empty(nsteps + 1)}
    hB = [[h * b for b in row] for row in _B]

    for k in range(nsteps + 1):
        f0, mem_in, s = dyn.evaluate(x, mem)
        sig["dfsys"][k] = s["dfsys"]
        sig["dfw"][k] = s["dfw"]
        if k % save_every == 0:
            for g in groups:
                out[g][k // save_every] = x[dyn.layout[g]]
        if k == nsteps:
            break
        f1, _, _ = dyn.evaluate(x + f0 * hB[0][0], mem)
        f2, _, _ = dyn.evaluate(x + (f0 * hB[1][0] + f1 * hB[1][1]), mem)
        x = x + (f0 * hB[2][0] + f1 * hB[2][1] + f2 * hB[2][2])
        mem = mem_in

    t = np.arange(nsave) * (h * save_every)
    return Result(t, out, np.arange(nsteps + 1) * h, sig)
