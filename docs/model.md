# The model

WAMCES models the continental European power system for studies of the frequency after a
disturbance. The network has 3809 buses and 7343 lines. Its sources are 618 synchronous machines,
300 grid-forming converters and 500 grid-following converters. A wide-area measurement system
of 1573 PMUs at the 400 kV buses reports to a central unit. The model and its data are by Musca,
Ippolito and Riva Sanseverino ([doi:10.3390/electricity7020028](https://doi.org/10.3390/electricity7020028),
[doi:10.5281/zenodo.18721887](https://doi.org/10.5281/zenodo.18721887)).

This page gives the equations of the Simulink model as the Python code implements them. Each
block of `wamces_model.slx` appears once in `Dynamics.evaluate` (`src/wamces/dynamics.py`). A
comment there names the block. The values of all gains and time constants are in [gains.md](gains.md). All quantities are per unit. The base of a machine or converter is its
own rating. The base of the network is 100 MVA.

## Initialisation

`Model.initialise` ports the script `wamces_script.m`.

1. The admittance matrix of the lines uses the series impedance `R + jX` of each line. The line
   charging `B` and the taps of the table are not used.
2. Each load becomes a constant admittance `(P - jQ) / V^2` at the voltage of the power flow.
   The source admittances of the machines and converters are added on the diagonal.
3. A machine that the power flow loads above 80 % of its rating gets the rating `S / 0.8`.
4. Kron reduction keeps the buses with a source (1418 buses). The matrix `Zred` maps the source
   currents to the source voltages. The matrix `Yaux` gives the voltages of the other buses.
5. The disturbance is a constant admittance `(dP / 100) / V^2` at the disturbance bus. It is
   part of the reduced matrices `Zredd` and `Yauxd` from the start. The simulation therefore
   starts with the disturbance in place.

The model uses the voltages of the 1030 PMU buses without a source. The voltages of the other
buses without a source are not computed.

## Network

The sources inject the currents `I` in the common xy frame of the network. The source voltages
are `V_inj = Zredd @ I_inj`. The voltages at the other PMU buses are `V_pmu = Yauxd_pmu @ V_inj`. A machine or converter turns a phasor `z` from the xy frame into its
own dq frame with `z * exp(-j delta)`, where `delta` is its angle.

## Synchronous machines

Each machine is a two-axis model with the transient voltages `Eq'` and `Ed'` (`sm_eq1`, `sm_ed1`).
With the terminal voltage `Vg` in the dq frame:

```
E      = Eq' + j (Ed' - (Xq' - Xd') iq_mem)
Idq    = Yg (E - Vg),   iq = Re Idq,   id = Im Idq
Pe     = Re(Vg conj(Idq))
dEq'   = (Efd + (Xd - Xd') id - Eq') / Td0'
dEd'   = (-(Xq - Xq') iq - Ed') / Tq0'
ddw    = (Pm - Pe - D dw) / M,   M = 2 H
ddelta = wn (dw - dfsys)
```

`iq_mem` is the output of a Memory block. It holds `iq` of the last major time step. The
machine injects the Norton current `Yg E` into the network. The angles follow the mean speed
deviation `dfsys` of all machines, so `delta` is an angle relative to that mean.

**Voltage regulator** (`avr_vt`, `avr_leadlag`, `avr_lag`). A lag `Tr` measures `|Vg|`. A
lead-lag `(1 + s TA) / (1 + s TB)` and the gain `KA` with the lag `TE` act on
`Vref + Vs - Vt`. The field voltage is that output plus `Efd0`. The model has no limits.

**Power system stabiliser** (`pss_*`). The speed deviation passes two washouts `TW1` and `TW2`.
The electrical power passes a washout `TW3`, the gain `KS2` and a lag `T7`. The sum
`KS1 (y_w + (KS3 - 1) y_p)` passes two lead-lags `(T1, T2)` and `(T3, T4)` and gives `Vs`.

**Governor and turbine** (`gov_lag`, `gov_leadlag`). The droop `(1 - w) / R` passes a lag `T1`
and a lead-lag `(T2, T3)`. The mechanical power is that output plus `Pm0`.

## Grid-forming converters

Each converter is a voltage source behind the filter admittance `Yf`. With the terminal
voltage `Vdq` in its own frame:

```
Vctrl   = Kvp (Vref - Vt) + xVctrl,   dxVctrl = Kvi (Vref - Vt),   dVt = (|Vdq| - Vt) / Tr
Pe      = Re(Vdq conj(Yf (Vctrl - Vdq)))
ddw     = Ki (Kw (dfw_d - dw) - Kf dw + P0 - Pe)
ddelta  = wn (dw - dfsys)
```

The converter injects `Yf Vctrl` in its own frame. `dfw_d` is the wide-area frequency `dfw`
after a latency `Td` (`gfm_lat1`, `gfm_lat2`). The latency is the second-order filter
`(6/Td^2 - (2/Td) s) / (s^2 + (4/Td) s + 6/Td^2)`. This is the Pade approximation of order
(1, 2) of a delay `Td`. The published script sets `Kw` to zero. The
wide-area channel then has no effect, but its states still move.

## Grid-following converters

Each converter has a phase-locked loop and a PI controller of its power.

```
dx1     = (900 / wn) Im Vdq
ddelta  = wn (x1 + (60 / wn) Im Vdq - dfsys)
id      = Kpp (P0 - P) + xP,   dxP = Kpi (P0 - P)
iq      = Kqp (Q0 - Q) + xQ,   dxQ = Kqi (Q0 - Q)
```

The converter injects `Yf Vdq_mem + id + j iq` in its own frame. Two Memory blocks hold
`Vdq_mem` and the measured power `P + jQ` of the last major time step. The measured power is
the product of the bus voltage in the xy frame and the current `id + j iq`, as in the block
`Product1`.

In the published power flow every grid-following converter has zero power. The controllers
then stay at rest. The PLL gains 900 and 60 are fixed numbers in the model. The columns `Kpllp`
and `Kplli` of the table hold the same values.

## PMUs and the central unit

Each PMU tracks the angle of its bus voltage with a PLL. The PI controller `(Kp / wn) + (Ki / wn) / s`
gives the speed deviation `dw` from `Im Vdq`. A lag of 0.04 s gives the measured frequency
`f = fn (1 + dw)` (`pmu_f`). The central unit passes each measured frequency through the latency
`Td` of its PMU. Its output `dfw` is the mean of the delayed frequencies divided by `fn`, minus one.

## Solver

The model is simulated with the fixed-step solver ode3 of Bogacki and Shampine. The step is
0.01 s. Simulink chooses this solver for the published settings. A Memory block keeps its output over a
whole step and gives the new value from the next major time step on. The details and the
comparison with Simulink are in [validation.md](validation.md).

## State groups

| Name | States | Simulink block |
|---|---|---|
| `sm_domega`, `sm_delta` | speed deviation, angle | Synchronous machines/mechanical |
| `sm_eq1`, `sm_ed1` | transient voltages | Synchronous machines/electrical |
| `avr_vt`, `avr_leadlag`, `avr_lag` | measured voltage, lead-lag, exciter | AVR |
| `gov_lag`, `gov_leadlag` | lag, lead-lag | Governor/Turbine |
| `pss_wash1`, `pss_wash2`, `pss_washp` | washouts of speed and power | PSS |
| `pss_lagp`, `pss_leadlag1`, `pss_leadlag2` | lag of power, lead-lags | PSS |
| `gfm_delta`, `gfm_domega` | angle, speed deviation | GFM generation/synchronization loop |
| `gfm_lat1`, `gfm_lat2` | latency of the wide-area frequency | GFM generation/synchronization loop/Latency |
| `gfm_vt`, `gfm_vctrl` | measured voltage, PI state | GFM generation/voltage controller |
| `gfl_pll1`, `gfl_delta` | PLL integrator, angle | GFL generation/synchronization loop |
| `gfl_ip`, `gfl_iq` | PI states of P and Q | GFL generation/P controller, Q controller |
| `pmu_theta`, `pmu_pi`, `pmu_f` | PLL angle, PI state, measured frequency | WAMC/PMUs |
| `cu_lat1`, `cu_lat2` | latency of the measured frequencies | WAMC/Central unit/Latency |

`wamces.dynamics.SIMULINK_BLOCKS` maps the full path of each Integrator block to these names.
