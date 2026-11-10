"""The block diagram of ``wamces_model.slx`` as one function of the state.

Every block of the Simulink model appears once in ``Dynamics.evaluate``. The comments name
the subsystem and the block. Each sum keeps the order of the inputs of its Simulink block.
Each gain is computed once from the tables, as Simulink evaluates its parameter.

The model has continuous states and three Memory blocks. A Memory block holds the value
of its input from the last major time step. The memories here hold the current q-axis
current of each machine, the dq voltage of each grid-following converter and the power
that its controllers measure.
"""

import numpy as np

from .model import FN, WN

# Order and names of the state groups in the flat state vector.
STATE_GROUPS = (
    "sm_domega", "sm_delta", "sm_eq1", "sm_ed1",
    "avr_vt", "avr_leadlag", "avr_lag",
    "gov_lag", "gov_leadlag",
    "pss_wash1", "pss_wash2", "pss_washp", "pss_lagp", "pss_leadlag1", "pss_leadlag2",
    "gfm_delta", "gfm_domega", "gfm_lat1", "gfm_lat2", "gfm_vt", "gfm_vctrl",
    "gfl_pll1", "gfl_delta", "gfl_ip", "gfl_iq",
    "pmu_theta", "pmu_pi", "pmu_f", "cu_lat1", "cu_lat2",
)
MEMORY_GROUPS = ("sm_iq", "gfl_vdq", "gfl_s")

# The Integrator block of each state group in wamces_model.slx.
_SG = "Sources/Synchronous generation/"
SIMULINK_BLOCKS = {
    _SG + "Synchronous machines/mechanical/xdomega": "sm_domega",
    _SG + "Synchronous machines/mechanical/xdelta": "sm_delta",
    _SG + "Synchronous machines/electrical/xEq1": "sm_eq1",
    _SG + "Synchronous machines/electrical/xEd1": "sm_ed1",
    _SG + "AVR/Lag (V)/int": "avr_vt",
    _SG + "AVR/Lead-lag/int": "avr_leadlag",
    _SG + "AVR/Lag/int": "avr_lag",
    _SG + "Governor//Turbine/Lag/int": "gov_lag",
    _SG + "Governor//Turbine/Lead-lag/int": "gov_leadlag",
    _SG + "PSS/Der-lag (dw) 1/PSS1": "pss_wash1",
    _SG + "PSS/Der-lag (dw) 2/PSS1": "pss_wash2",
    _SG + "PSS/Der-lag (P) 1/PSS1": "pss_washp",
    _SG + "PSS/Lag (P)/int": "pss_lagp",
    _SG + "PSS/Lead-lag 1/PSS2": "pss_leadlag1",
    _SG + "PSS/Lead-lag 2/PSS3": "pss_leadlag2",
    "Sources/GFM generation/synchronization loop/xdelta": "gfm_delta",
    "Sources/GFM generation/synchronization loop/xdomega": "gfm_domega",
    "Sources/GFM generation/synchronization loop/Latency/Integrator1": "gfm_lat1",
    "Sources/GFM generation/synchronization loop/Latency/Integrator2": "gfm_lat2",
    "Sources/GFM generation/voltage controller/Vtlag": "gfm_vt",
    "Sources/GFM generation/voltage controller/xVctrl": "gfm_vctrl",
    "Sources/GFL generation/synchronization loop/x1PLL": "gfl_pll1",
    "Sources/GFL generation/synchronization loop/x2PLL": "gfl_delta",
    "Sources/GFL generation/P controller/xPerr": "gfl_ip",
    "Sources/GFL generation/Q controller/xPerr": "gfl_iq",
    "WAMC/PMUs/Integrator": "pmu_theta",
    "WAMC/PMUs/PID Controller/Integrator/Continuous/Integrator": "pmu_pi",
    "WAMC/PMUs/LPF/int": "pmu_f",
    "WAMC/Central unit/Latency/Integrator1": "cu_lat1",
    "WAMC/Central unit/Latency/Integrator2": "cu_lat2",
}


def cplx(re, im):
    """The complex array re + 1i*im, built without rounding."""
    out = np.empty(np.shape(re), dtype=complex)
    out.real = re
    out.imag = im
    return out


def seqsum(v):
    """Sum of the elements in index order, as the Simulink block Sum of Elements."""
    return np.cumsum(v)[-1]


class Dynamics:
    """Right-hand side of the WAMCES model.

    ``evaluate(x, mem)`` returns the time derivative of the state, the new inputs of the
    Memory blocks and a dict of signals. The state ``x`` is one flat array. The method
    ``unpack`` splits it into the named groups of ``STATE_GROUPS``.
    """

    def __init__(self, model):
        self.model = model
        sizes = {k: len(v) for k, v in model.x0.items()}
        self.layout = {}
        start = 0
        for name in STATE_GROUPS:
            self.layout[name] = slice(start, start + sizes[name])
            start += sizes[name]
        self.nx = start
        self._gains(model)

    # ------------------------------------------------------------------ helpers
    def pack(self, groups):
        x = np.empty(self.nx)
        for name, sl in self.layout.items():
            x[sl] = groups[name]
        return x

    def unpack(self, x):
        return {name: x[sl] for name, sl in self.layout.items()}

    def _gains(self, m):
        sm, avr, tur, pss, gfm, gfl, pmu = m.sm, m.avr, m.tur, m.pss, m.gfm, m.gfl, m.pmu
        g = self
        # Positions of the sources inside the reduced network.
        pos = {b: k for k, b in enumerate(m.index_inj)}
        g.pos_sg = np.array([pos[b] for b in m.index_sg])
        g.pos_gfm = np.array([pos[b] for b in m.index_gfm])
        g.pos_gfl = np.array([pos[b] for b in m.index_gfl])
        g.Zredd = np.ascontiguousarray(m.Zredd)
        g.Yauxd_pmu = np.ascontiguousarray(m.Yauxd_pmu)
        g.nbus = m.nbus

        # Synchronous machines: electrical and mechanical part.
        g.sm_inv_td01 = 1.0 / sm.Td01
        g.sm_inv_tq01 = 1.0 / sm.Tq01
        g.sm_xd_xd1 = sm.Xd - sm.Xd1
        g.sm_xq_xq1 = sm.Xq - sm.Xq1
        g.sm_xq1_xd1 = sm.Xq1 - sm.Xd1
        g.sm_yg = sm.Yg
        g.sm_inv_m = 1.0 / sm.M
        g.sm_d = sm.D
        g.sm_sr2sb = sm.Sr2Sb
        g.sm_vref = sm.Vref
        g.sm_efd0 = sm.Efd0
        g.sm_pm0 = sm.Pm0
        # AVR
        g.avr_inv_tr = 1.0 / avr.Tr
        g.avr_inv_tb = 1.0 / avr.TB
        g.avr_tab = (avr.TA - avr.TB) / avr.TB
        g.avr_ka = avr.KA
        g.avr_inv_te = 1.0 / avr.TE
        # Governor and turbine
        g.tur_inv_r = 1.0 / tur.R
        g.tur_inv_t1 = 1.0 / tur.T1
        g.tur_inv_t3 = 1.0 / tur.T3
        g.tur_t23 = (tur.T2 - tur.T3) / tur.T3
        # PSS
        g.pss_inv_tw1 = 1.0 / pss.TW1
        g.pss_inv_tw2 = 1.0 / pss.TW2
        g.pss_inv_tw3 = 1.0 / pss.TW3
        g.pss_inv_t7 = 1.0 / pss.T7
        g.pss_ks1, g.pss_ks2, g.pss_ks3 = pss.KS1, pss.KS2, pss.KS3
        g.pss_inv_t2 = 1.0 / pss.T2
        g.pss_t12 = (pss.T1 - pss.T2) / pss.T2
        g.pss_inv_t4 = 1.0 / pss.T4
        g.pss_t34 = (pss.T3 - pss.T4) / pss.T4

        # Grid-forming converters
        g.gfm_yf = gfm.Yf
        g.gfm_sr2sb = gfm.Sr2Sb
        g.gfm_vref = gfm.Vref
        g.gfm_kvp, g.gfm_kvi = gfm.Kvp, gfm.Kvi
        g.gfm_inv_tr = 1.0 / gfm.Tr
        g.gfm_kf, g.gfm_ki, g.gfm_kw = gfm.Kf, gfm.Ki, gfm.Kw
        g.gfm_p0 = gfm.P0
        g.gfm_lat_g1 = 4.0 / gfm.Td
        g.gfm_lat_g2 = 6.0 / gfm.Td**2
        g.gfm_lat_g3 = 2.0 / gfm.Td

        # Grid-following converters. The PLL gains are fixed numbers in the model.
        g.gfl_yf = gfl.Yf
        g.gfl_sr2sb = gfl.Sr2Sb
        g.gfl_p0, g.gfl_q0 = gfl.P0, gfl.Q0
        g.gfl_kpp, g.gfl_kpi = gfl.Kpp, gfl.Kpi
        g.gfl_kqp, g.gfl_kqi = gfl.Kqp, gfl.Kqi
        g.gfl_pll_ki = 900 / WN
        g.gfl_pll_kp = 60 / WN

        # PMUs and the central unit
        g.pmu_idx = m.index_pmu
        g.pmu_p = pmu.Kp / WN
        g.pmu_i = pmu.Ki / WN
        g.pmu_inv_tf = 1.0 / 0.04
        g.cu_g1 = 4.0 / pmu.Td
        g.cu_g2 = 6.0 / pmu.Td**2
        g.cu_g3 = 2.0 / pmu.Td
        g.cu_n = len(m.index_pmu)
        g.index_inj = m.index_inj
        g.index_gfm = m.index_gfm
        g.index_gfl = m.index_gfl
        g.index_pmu_bs = m.index_pmu_bs

    # ------------------------------------------------------------------ the model
    def evaluate(self, x, mem):
        """Derivatives, Memory inputs and signals at the state ``x``.

        ``mem`` holds the outputs of the three Memory blocks.
        """
        g = self
        s = self.unpack(x)
        dx = np.empty(self.nx)
        d = {name: dx[sl] for name, sl in self.layout.items()}

        # ---- Synchronous generation/Synchronous machines/mechanical: speeds
        domega = s["sm_domega"]
        dfsys = seqsum(domega) / len(domega)  # Sum of Elements, Divide by Width
        w = domega + 1.0  # Add (xdomega, Constant)

        # ---- WAMC/Central unit: the latency of the measured frequencies
        c1, c2 = s["cu_lat1"], s["cu_lat2"]
        y_cu = g.cu_g2 * c2 - g.cu_g3 * c1  # Latency/Sum1 (-+)
        dfw = (seqsum(y_cu) / g.cu_n) * (1 / FN) - 1.0  # Divide, Gain2, Add2 (+-)

        # ---- Sources: the current injections of all sources
        # Synchronous machines/electrical
        delta = s["sm_delta"]
        iq_mem = mem["sm_iq"]
        E = cplx(s["sm_eq1"], s["sm_ed1"] - g.sm_xq1_xd1 * iq_mem)  # Add3 (+-), Real-Imag to Complex1
        Ig = g.sm_yg * E  # Gain9
        I_sm = g.sm_sr2sb * (Ig * np.exp(1j * delta))  # dq to xy, Sr -> Sb

        # GFM generation
        gd = s["gfm_delta"]
        a1 = g.gfm_vref - s["gfm_vt"]  # voltage controller/Add1 (+-)
        vctrl = g.gfm_kvp * a1 + s["gfm_vctrl"]  # voltage controller/Add
        E_gfm = cplx(vctrl, np.zeros_like(vctrl))  # Real-Imag to Complex
        I_gfm = g.gfm_sr2sb * ((g.gfm_yf * E_gfm) * np.exp(1j * gd))  # Gain2, dq to xy, Sr -> Sb

        # GFL generation
        ld = s["gfl_delta"]
        S_mem = mem["gfl_s"]
        eP = g.gfl_p0 - S_mem.real  # P controller/Sum2 (+-)
        eQ = g.gfl_q0 - S_mem.imag  # Q controller/Sum2 (+-)
        i_d = g.gfl_kpp * eP + s["gfl_ip"]  # P controller/Sum7
        i_q = g.gfl_kqp * eQ + s["gfl_iq"]  # Q controller/Sum7
        I_ctrl = cplx(i_d, i_q)
        I_gfl = g.gfl_sr2sb * ((g.gfl_yf * mem["gfl_vdq"] + I_ctrl) * np.exp(1j * ld))  # Add2, dq to xy

        # ---- Network: the source voltages, then the voltages at the other PMU buses
        I_inj = np.empty(len(g.index_inj), dtype=complex)
        I_inj[g.pos_sg] = I_sm
        I_inj[g.pos_gfm] = I_gfm
        I_inj[g.pos_gfl] = I_gfl
        V_inj = g.Zredd @ I_inj  # Zred
        Vbus = np.zeros(g.nbus, dtype=complex)  # Constant5, Real-Imag to Complex
        Vbus[g.index_inj] = V_inj  # Assign (src)
        Vbus[g.index_pmu_bs] = g.Yauxd_pmu @ V_inj  # Yaux, Assign (pmu)

        # ---- Synchronous machines: electrical part
        Vg = V_inj[g.pos_sg] * np.exp(-1j * delta)  # xy to dq
        Idq = g.sm_yg * (E - Vg)  # Add7 (+-), Gain10
        iq, i_dsm = Idq.real, Idq.imag
        pe = Vg.real * iq + Vg.imag * i_dsm  # Product with Conj, real part
        vt = np.abs(Vg)
        # PSS, then the AVR
        vs = self._pss(s, d, w, pe)
        u_avr = (g.sm_vref + vs) - s["avr_vt"]  # AVR/Sum (++-)
        d["avr_vt"][:] = g.avr_inv_tr * (vt - s["avr_vt"])  # Lag (V)
        e_ll = u_avr - s["avr_leadlag"]  # Lead-lag/Sum5
        d["avr_leadlag"][:] = g.avr_inv_tb * e_ll
        y_ll = u_avr + g.avr_tab * e_ll  # Lead-lag/Sum3
        d["avr_lag"][:] = g.avr_inv_te * (g.avr_ka * y_ll - s["avr_lag"])  # Gain5, Lag
        efd = s["avr_lag"] + g.sm_efd0  # AVR/Sum1
        # Governor and turbine
        u_gov = g.tur_inv_r * (1.0 - w)  # Sum (|-+), 1/R
        d["gov_lag"][:] = g.tur_inv_t1 * (u_gov - s["gov_lag"])  # Lag
        e_gov = s["gov_lag"] - s["gov_leadlag"]  # Lead-lag/Sum5
        d["gov_leadlag"][:] = g.tur_inv_t3 * e_gov
        pm = (s["gov_lag"] + g.tur_t23 * e_gov) + g.sm_pm0  # Lead-lag/Sum3, Sum1
        # electrical states
        d["sm_eq1"][:] = g.sm_inv_td01 * ((efd + g.sm_xd_xd1 * i_dsm) - s["sm_eq1"])  # Add (++-)
        d["sm_ed1"][:] = g.sm_inv_tq01 * (-(g.sm_xq_xq1 * iq) - s["sm_ed1"])  # Add1 (--)
        # mechanical states
        d["sm_domega"][:] = g.sm_inv_m * ((pm - pe) - g.sm_d * domega)  # Sub (+--), Gain3
        d["sm_delta"][:] = WN * (domega - dfsys)  # Add1 (-+), Gain

        # ---- GFM generation
        Vdq = Vbus[g.index_gfm] * np.exp(-1j * gd)  # xy to dq
        d["gfm_vt"][:] = g.gfm_inv_tr * (np.abs(Vdq) - s["gfm_vt"])  # voltage controller/Sum2
        d["gfm_vctrl"][:] = g.gfm_kvi * a1  # voltage controller/Gain1
        I_f = g.gfm_yf * (E_gfm - Vdq)  # Add2 (+-), Gain3
        pe_gfm = Vdq.real * I_f.real + Vdq.imag * I_f.imag  # Product with Conj, real part
        l1, l2 = s["gfm_lat1"], s["gfm_lat2"]
        d["gfm_lat1"][:] = (dfw - g.gfm_lat_g2 * l2) - g.gfm_lat_g1 * l1  # Latency/Sum (+--)
        d["gfm_lat2"][:] = l1
        y_lat = g.gfm_lat_g2 * l2 - g.gfm_lat_g3 * l1  # Latency/Sum1 (-+)
        dom = s["gfm_domega"]
        sub = ((g.gfm_kw * (y_lat - dom) - g.gfm_kf * dom) + g.gfm_p0) - pe_gfm  # Sum2, Sub (+-+-)
        d["gfm_domega"][:] = g.gfm_ki * sub  # Ki
        d["gfm_delta"][:] = WN * (dom - dfsys)  # Add1 (-+), Gain

        # ---- GFL generation
        V_gfl = Vbus[g.index_gfl]
        Vdq_l = V_gfl * np.exp(-1j * ld)  # xy to dq
        vi = Vdq_l.imag
        d["gfl_pll1"][:] = g.gfl_pll_ki * vi  # Gain5
        d["gfl_delta"][:] = WN * ((s["gfl_pll1"] + g.gfl_pll_kp * vi) - dfsys)  # Add4, sum_omega
        d["gfl_ip"][:] = g.gfl_kpi * eP  # P controller/Gain
        d["gfl_iq"][:] = g.gfl_kqi * eQ  # Q controller/Gain

        # ---- WAMC/PMUs
        th = s["pmu_theta"]
        Vp = Vbus[g.pmu_idx] * np.exp(-1j * th)  # xy to dq
        ep = Vp.imag
        dw = g.pmu_p * ep + s["pmu_pi"]  # PID Controller, parallel form
        d["pmu_pi"][:] = g.pmu_i * ep
        d["pmu_theta"][:] = WN * (dw - dfsys)  # Add1 (-+), Gain1
        d["pmu_f"][:] = g.pmu_inv_tf * (FN * (dw + 1.0) - s["pmu_f"])  # Add, Gain, LPF
        # Central unit/Latency
        d["cu_lat1"][:] = (s["pmu_f"] - g.cu_g2 * c2) - g.cu_g1 * c1  # Latency/Sum (+--)
        d["cu_lat2"][:] = c1

        mem_in = {"sm_iq": iq, "gfl_vdq": Vdq_l, "gfl_s": V_gfl * I_ctrl}
        sig = {"dfsys": dfsys, "dfw": dfw, "Vbus": Vbus}
        return dx, mem_in, sig

    def _pss(self, s, d, w, pe):
        """Power system stabilisers. Returns the stabilising voltage Vs."""
        g = self
        u = w - 1.0  # Sum6 (|+-)
        y1 = u - s["pss_wash1"]  # Der-lag (dw) 1
        d["pss_wash1"][:] = g.pss_inv_tw1 * y1
        y2 = y1 - s["pss_wash2"]  # Der-lag (dw) 2
        d["pss_wash2"][:] = g.pss_inv_tw2 * y2
        yp = pe - s["pss_washp"]  # Der-lag (P) 1
        d["pss_washp"][:] = g.pss_inv_tw3 * yp
        d["pss_lagp"][:] = g.pss_inv_t7 * (g.pss_ks2 * yp - s["pss_lagp"])  # Gain1, Lag (P)
        lp = s["pss_lagp"]
        u1 = g.pss_ks1 * ((y2 + g.pss_ks3 * lp) - lp)  # Add, Add1 (+-), Gain
        e1 = u1 - s["pss_leadlag1"]  # Lead-lag 1
        d["pss_leadlag1"][:] = g.pss_inv_t2 * e1
        z1 = u1 + g.pss_t12 * e1
        e2 = z1 - s["pss_leadlag2"]  # Lead-lag 2
        d["pss_leadlag2"][:] = g.pss_inv_t4 * e2
        return z1 + g.pss_t34 * e2
