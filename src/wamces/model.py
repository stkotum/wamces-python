"""The WAMCES case: parameters, network reduction and initial conditions.

This module ports the published initialisation script ``wamces_script.m`` line by line.
Each step keeps the order of operations of the script. The rounding of the results
therefore stays close to MATLAB's. Indices are 0-based here and 1-based in the script.
"""

import numpy as np

from .data import Table, load_tables

SB = 100.0  # system base power in MVA
FN = 50.0  # nominal frequency in Hz
WN = 2 * np.pi * FN  # nominal angular frequency in rad/s

# The disturbance of the published script: a load increase of 1000 MW at bus 2178.
DISTURBANCE_BUS = 2178
DISTURBANCE_MW = 1000.0

# A machine that the power flow loads above this fraction of its rating gets a larger rating.
MAX_LOADING = 0.8


def admittance_matrix(nbus, lin):
    """Bus admittance matrix of the lines as a dict that maps (i, j) to y_ij.

    This is ``build_admittance_matrix`` of the script. The entries are summed in the
    order of the line table, as in the script. Line charging and taps are not used.
    """
    Y = {}
    for i, j, r, x in zip(lin.ID1.astype(int) - 1, lin.ID2.astype(int) - 1, lin.R, lin.X):
        z = complex(r, x)
        for key, sign in (((i, j), -1), ((j, i), -1), ((i, i), 1), ((j, j), 1)):
            y = 1 / z
            old = Y.get(key, 0j)
            Y[key] = old - y if sign < 0 else old + y
    return Y


def dense_block(Y, rows, cols):
    """The block Y(rows, cols) of an admittance dict as a dense complex array."""
    pos_r = {b: k for k, b in enumerate(rows)}
    pos_c = {b: k for k, b in enumerate(cols)}
    out = np.zeros((len(rows), len(cols)), dtype=complex)
    for (i, j), y in Y.items():
        if i in pos_r and j in pos_c:
            out[pos_r[i], pos_c[j]] = y
    return out


def reduce_network(Y_aug, index_bs, index_inj, Yaux=None):
    """Kron reduction of the network onto the buses with sources.

    Returns ``Yaux`` and ``Yred`` of the script. The voltages of the buses without a
    source follow from the source voltages by ``V_bs = Yaux @ V_inj``. A given ``Yaux``
    is used as it is.
    """
    if Yaux is None:
        A = dense_block(Y_aug, index_bs, index_bs)
        B = dense_block(Y_aug, index_bs, index_inj)
        Yaux = np.linalg.solve(-A, B)
    C = dense_block(Y_aug, index_inj, index_bs)
    D = dense_block(Y_aug, index_inj, index_inj)
    return Yaux, D + C @ Yaux


def _phasor(Vm, Va):
    return Vm * (np.cos(Va) + 1j * np.sin(Va))


def _rotate(z, delta):
    """z times exp(-1i*delta), the rotation from the xy frame into the dq frame."""
    return z * np.exp(-1j * delta)


class Model:
    """The WAMCES model of the continental European power system.

    The constructor runs the initialisation of the published script. The disturbance is
    a load step at one bus. It acts from the start of the simulation, as in the
    published model.

    The tables of the machines, controllers and converters are attributes. A change of
    a table column followed by ``initialise`` gives a new operating point. The rating of
    a machine in the data is ``sm.Sr_data``. The column ``sm.Sr`` holds the rating after
    the rule for overloaded machines, as ``Tsm.Sr`` in the script. A change of a
    derived value such as ``gfl.P0`` or ``gfm.Kw`` after ``initialise`` changes the
    dynamics only.

    Parameters
    ----------
    tables : dict of Table, optional
        The output of ``load_tables``. The default reads the files of the package.
    disturbance_bus : int
        The bus of the load step, 1-based as in the data.
    disturbance_mw : float
        The size of the load step in MW. Zero gives the undisturbed operating point.
    kw : float or array
        The gain of the wide-area channel of the grid-forming converters. The published
        script sets it to zero.
    """

    def __init__(self, tables=None, disturbance_bus=DISTURBANCE_BUS,
                 disturbance_mw=DISTURBANCE_MW, kw=0.0):
        t = load_tables() if tables is None else tables
        self.bus = t["bus"].copy()
        self.lin = t["lin"].copy()
        self.sm = t["sm"].copy()
        self.avr = t["avr"].copy()
        self.pss = t["pss"].copy()
        self.tur = t["tur"].copy()
        self.gfm = t["gfm"].copy()
        self.gfl = t["gfl"].copy()
        self.pmu = t["pmu"].copy()
        self.pf = t["pf"].copy()
        # The script raises the rating of overloaded machines in Tsm.Sr. The rating of
        # the data is kept apart, so that each call of initialise starts from it.
        if "Sr_data" not in self.sm:
            self.sm.Sr_data = self.sm.Sr.copy()
        self.disturbance_bus = int(disturbance_bus)
        self.disturbance_mw = float(disturbance_mw)
        self.kw = kw
        self.initialise()

    # ------------------------------------------------------------------ initialisation
    def initialise(self):
        """Run the initialisation of the published script on the current tables."""
        self._sizes_and_indices()
        self._machines()
        self._converters()
        self._network()
        self._initial_state()

    def _sizes_and_indices(self):
        bus = self.bus
        self.nbus = bus.nrows
        self.nsm = self.sm.nrows
        self.nlin = self.lin.nrows
        self.index_sg = self.sm.BusID.astype(int) - 1
        self.index_gfm = self.gfm.BusID.astype(int) - 1
        self.index_gfl = self.gfl.BusID.astype(int) - 1
        index_bs = np.flatnonzero(bus.Type == 3)
        self.index_400 = np.flatnonzero(bus.Vr > 230)
        self.ngfm = len(self.index_gfm)
        self.ngfl = len(self.index_gfl)
        self.ninj = self.nsm + self.ngfm + self.ngfl
        self.index_ibr = np.union1d(self.index_gfm, self.index_gfl)
        self.index_inj = np.union1d(self.index_sg, self.index_ibr)
        self.index_bs = np.setdiff1d(index_bs, self.index_ibr)
        self.index_pmu = self.index_400
        self.index_pmu_bs = np.setdiff1d(self.index_pmu, self.index_inj)
        if len(self.index_inj) != self.ninj:
            raise ValueError("two sources share a bus")
        if len(self.index_pmu) != self.pmu.nrows:
            raise ValueError("the PMU table does not match the 400 kV buses")

        # The columns of the power-flow solution that the script uses.
        pf = self.pf
        self.Vm = pf.Vmag
        self.Va = pf.Vang

    def _machines(self):
        sm, pf = self.sm, self.pf
        ig = self.index_sg
        sm.P = pf.Pgen[ig]
        sm.Q = pf.Qgen[ig]
        sm.Sg = np.sqrt(sm.P**2 + sm.Q**2) * SB
        over = sm.Sg > sm.Sr_data * MAX_LOADING
        self.overloaded = np.flatnonzero(over)
        sm.Sr = np.where(over, sm.Sg / MAX_LOADING, sm.Sr_data)
        sm.Sr2Sb = sm.Sr / SB
        sm.M = 2 * sm.H
        sm.Yg = 1 / (sm.Ra + 1j * sm.Xd1)

        sm.Vm = pf.Vmag[ig]
        sm.Va = pf.Vang[ig]
        sm.P = pf.Pgen[ig] / sm.Sr2Sb
        sm.Q = pf.Qgen[ig] / sm.Sr2Sb
        sm.Vxy = _phasor(sm.Vm, sm.Va)
        sm.Ixy = np.conj((sm.P + 1j * sm.Q) / sm.Vxy)
        sm.Eq0 = sm.Vxy + (sm.Ra + 1j * sm.Xq) * sm.Ixy
        sm.delta0 = np.angle(sm.Eq0)
        sm.Iqd = _rotate(sm.Ixy, sm.delta0)
        sm.Iq0 = sm.Iqd.real.copy()
        sm.Id0 = sm.Iqd.imag.copy()
        sm.Vqd = _rotate(sm.Vxy, sm.delta0)
        sm.Vq0 = sm.Vqd.real.copy()
        sm.Vd0 = sm.Vqd.imag.copy()
        sm.Efd0 = np.abs(sm.Eq0) - (sm.Xd - sm.Xq) * sm.Id0
        sm.Eq_dash0 = sm.Efd0 + (sm.Xd - sm.Xd1) * sm.Id0
        sm.Ed_dash0 = -(sm.Xq - sm.Xq1) * sm.Iq0
        sm.Pe0 = sm.Vq0 * sm.Iq0 + sm.Vd0 * sm.Id0
        sm.Eqd_dash = sm.Eq_dash0 + 1j * sm.Ed_dash0
        sm.Eqd_internal = sm.Eqd_dash - 1j * (sm.Xq1 - sm.Xd1) * sm.Iq0
        sm.Igxy = sm.Yg * sm.Vxy + sm.Ixy
        sm.Vref = sm.Vm.copy()
        sm.Pm0 = sm.Pe0.copy()

    def _converters(self):
        pf = self.pf
        self.gfm.Kw = np.zeros(self.ngfm) + self.kw
        for c, idx in ((self.gfm, self.index_gfm), (self.gfl, self.index_gfl)):
            c.Sr2Sb = c.Sr / SB
            c.Vm = pf.Vmag[idx]
            c.Va = pf.Vang[idx]
            c.P = pf.Pgen[idx] / c.Sr2Sb
            c.Q = pf.Qgen[idx] / c.Sr2Sb
            c.Vxy = _phasor(c.Vm, c.Va)
            c.Ixy = np.conj((c.P + 1j * c.Q) / c.Vxy)
            c.Exy_internal = c.Vxy + c.Ixy / c.Yf
            c.delta0 = np.angle(c.Exy_internal)
            c.Idq = _rotate(c.Ixy, c.delta0)
            c.Vdq = _rotate(c.Vxy, c.delta0)
            c.Edq_internal = _rotate(c.Exy_internal, c.delta0)
            c.Igxy = c.Yf * c.Vxy + c.Ixy
            c.Vref = c.Vm.copy()
            c.P0 = c.P.copy()
            c.Q0 = c.Q.copy()

    def _network(self):
        nbus, pf = self.nbus, self.pf
        Y = admittance_matrix(nbus, self.lin)
        self.Y = Y
        # The shunts of the loads and of the sources lie on the diagonal. They are
        # added in the order of the script: Y + Yload + Ysm + Ygfm + Ygfl.
        diag = {
            "load": (pf.Ploa - 1j * pf.Qloa) / self.Vm**2,
            "sm": np.zeros(nbus, dtype=complex),
            "gfm": np.zeros(nbus, dtype=complex),
            "gfl": np.zeros(nbus, dtype=complex),
        }
        diag["sm"][self.index_sg] = self.sm.Yg * self.sm.Sr2Sb
        diag["gfm"][self.index_gfm] = self.gfm.Yf * self.gfm.Sr2Sb
        diag["gfl"][self.index_gfl] = self.gfl.Yf * self.gfl.Sr2Sb
        Y_aug = dict(Y)
        for i in range(nbus):
            v = Y_aug.get((i, i), 0j)
            for key in ("load", "sm", "gfm", "gfl"):
                v = v + diag[key][i]
            Y_aug[(i, i)] = v
        self.Y_aug = Y_aug

        self.Yaux, self.Yred = reduce_network(Y_aug, self.index_bs, self.index_inj)
        self._Zred = None

        # The load step enters as a constant admittance at the disturbance bus.
        idis = self.disturbance_bus - 1
        Y_augd = dict(Y_aug)
        Y_augd[(idis, idis)] = Y_augd[(idis, idis)] + (self.disturbance_mw / 100) / (self.Vm[idis] ** 2)
        self.Y_augd = Y_augd
        # A step at a bus with a source leaves the blocks of the other buses unchanged.
        # The script computes Yauxd anyway and gets the same matrix as Yaux.
        same = idis in set(self.index_inj.tolist())
        self.Yauxd, self.Yredd = reduce_network(Y_augd, self.index_bs, self.index_inj,
                                                self.Yaux if same else None)
        self.Zredd = np.linalg.solve(self.Yredd, np.eye(self.ninj))

        pos = {b: k for k, b in enumerate(self.index_bs)}
        self.index_pmu_bs_yaux = np.array([pos[b] for b in self.index_pmu_bs], dtype=int)
        self.Yaux_pmu = self.Yaux[self.index_pmu_bs_yaux, :]
        self.Yauxd_pmu = self.Yauxd[self.index_pmu_bs_yaux, :]

    @property
    def Zred(self):
        """The inverse of ``Yred``, the reduced network without the disturbance.

        The simulation does not need it. It is computed on first use.
        """
        if self._Zred is None:
            self._Zred = np.linalg.solve(self.Yred, np.eye(self.ninj))
        return self._Zred

    def _initial_state(self):
        """Initial conditions of all integrators and memories of the Simulink model."""
        sm, gfm, gfl, pmu = self.sm, self.gfm, self.gfl, self.pmu
        n, ng, nl, npmu = self.nsm, self.ngfm, self.ngfl, len(self.index_pmu)
        z = np.zeros
        self.x0 = {
            # synchronous machines
            "sm_domega": z(n),
            "sm_delta": sm.delta0.copy(),
            "sm_eq1": sm.Eq_dash0.copy(),
            "sm_ed1": sm.Ed_dash0.copy(),
            "avr_vt": sm.Vm.copy(),
            "avr_leadlag": z(n),
            "avr_lag": z(n),
            "gov_lag": z(n),
            "gov_leadlag": z(n),
            "pss_wash1": z(n),
            "pss_wash2": z(n),
            "pss_washp": sm.Pe0.copy(),
            "pss_lagp": z(n),
            "pss_leadlag1": z(n),
            "pss_leadlag2": z(n),
            # grid-forming converters
            "gfm_delta": gfm.delta0.copy(),
            "gfm_domega": z(ng),
            "gfm_lat1": z(ng),
            "gfm_lat2": z(ng),
            "gfm_vt": gfm.Vm.copy(),
            "gfm_vctrl": np.abs(gfm.Edq_internal),
            # grid-following converters
            "gfl_pll1": z(nl),
            "gfl_delta": gfl.delta0.copy(),
            "gfl_ip": gfl.Idq.real.copy(),
            "gfl_iq": gfl.Idq.imag.copy(),
            # phasor measurement units and the central unit
            "pmu_theta": self.Va[self.index_pmu].copy(),
            "pmu_pi": z(npmu),
            "pmu_f": FN + z(npmu),
            "cu_lat1": z(npmu),
            "cu_lat2": FN / (6.0 / pmu.Td**2),
        }
        self.m0 = {
            "sm_iq": sm.Iq0.copy(),
            "gfl_vdq": gfl.Vdq.copy(),
            "gfl_s": gfl.P + 1j * gfl.Q,
        }

    # ------------------------------------------------------------------ simulation
    def simulate(self, t_end=50.0, h=0.01, save_every=10, save=None):
        """Simulate the model. See ``wamces.simulate.simulate``."""
        from .simulate import simulate

        return simulate(self, t_end=t_end, h=h, save_every=save_every, save=save)
