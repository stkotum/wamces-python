# Validation

The port was compared with the published Simulink model in MATLAB R2025b. The comparison covers
the data, the initialisation and every continuous state of the simulation at every step. All
numbers on this page come from the runs of 10 November 2026.

## Data

The CSV files of the package equal the tables of the published .mat files exactly. The CSV
files of the Zenodo record differ from them in one column. They round `Kf` of `data_gfm.csv` to
15 digits, which changes it by 5e-16 relative. `tools/make_data.py` rebuilds the files of the
package byte for byte from the output of `tools/export_tables.m`.

## Initialisation

The Python initialisation was compared with the workspace of the published script.

| Quantity | Agreement |
|---|---|
| all index vectors (sources, buses without a source, PMUs) | equal |
| ratings, per-unit powers, source admittances | equal |
| machine initial values (`delta0`, `Iq0`, `Efd0`, `Eq_dash0`, `Pe0` and others) | within 6e-16 relative |
| converter initial values (`delta0`, `Vdq`, `Edq_internal`) | within 4e-16 relative |
| augmented admittance matrix at the source buses | within 6e-19 relative |
| reduced matrices `Zredd` and `Yauxd_pmu` | within 6e-15 relative |

The differences of the reduced matrices come from the LU factorisation. MATLAB and numpy use
the same algorithm with different libraries. The small differences in the initial values come
from the complex exponential and the complex division.

## Simulation

`tools/make_reference.m` simulated two cases with every state logged at every step.

1. **Published**: the published script and model, 50 s, 5001 steps.
2. **Variant**: 10 s, 1001 steps. The wide-area gain is `Kw = 20`. The load step of 1000 MW
   is at bus 1. The power setpoints of the grid-following converters are raised by 0.1 pu.

The variant drives three parts of the model that rest in the published case. The wide-area gain
is zero in the published script. Bus 1 has no source. A step there changes `Yauxd` as well as
`Zredd`. The grid-following converters have zero power in the published power flow. Their power
controllers act only after a change of the setpoints.

Both cases have 20 935 continuous states. The table gives the largest absolute error over all
states of a group and all steps. The relative error divides it by the largest absolute value of
the group in Simulink.

| Group | States | Published: error | Published: relative | Variant: error | Variant: relative |
|---|---|---|---|---|---|
| `sm_domega` | 618 | 3.5e-16 | 3.0e-13 | 2.6e-16 | 5.5e-13 |
| `sm_delta` | 618 | 2.7e-13 | 1.3e-13 | 1.1e-13 | 6.5e-14 |
| `sm_eq1` | 618 | 2.3e-13 | 1.7e-13 | 9.3e-14 | 7.3e-14 |
| `sm_ed1` | 618 | 1.1e-13 | 1.7e-13 | 5.2e-14 | 8.3e-14 |
| `avr_vt` | 618 | 1.1e-14 | 1.1e-14 | 1.3e-14 | 1.3e-14 |
| `avr_leadlag` | 618 | 1.7e-15 | 5.8e-13 | 1.3e-15 | 4.5e-12 |
| `avr_lag` | 618 | 1.9e-12 | 5.2e-13 | 1.9e-12 | 1.9e-12 |
| `gov_lag` | 618 | 7.2e-16 | 2.3e-13 | 4.3e-16 | 3.0e-13 |
| `gov_leadlag` | 618 | 2.6e-16 | 1.2e-13 | 1.7e-16 | 2.9e-13 |
| `pss_wash1` | 618 | 1.7e-16 | 1.8e-13 | 9.5e-17 | 2.3e-13 |
| `pss_wash2` | 618 | 8.5e-17 | 2.5e-13 | 5.1e-17 | 4.4e-13 |
| `pss_washp` | 618 | 1.4e-15 | 1.8e-15 | 1.2e-15 | 1.5e-15 |
| `pss_lagp` | 618 | 1.8e-16 | 2.6e-13 | 1.3e-16 | 5.5e-13 |
| `pss_leadlag1` | 618 | 3.1e-15 | 3.0e-13 | 2.4e-15 | 6.2e-13 |
| `pss_leadlag2` | 618 | 1.1e-14 | 3.3e-13 | 1.0e-14 | 1.1e-12 |
| `gfm_delta` | 300 | 5.4e-14 | 2.6e-14 | 3.7e-14 | 2.3e-14 |
| `gfm_domega` | 300 | 3.7e-16 | 2.7e-13 | 1.5e-16 | 3.2e-13 |
| `gfm_lat1` | 300 | 3.6e-18 | 3.6e-12 | 3.2e-18 | 8.0e-12 |
| `gfm_lat2` | 300 | 2.5e-19 | 1.7e-13 | 2.2e-19 | 2.9e-13 |
| `gfm_vt` | 300 | 1.8e-14 | 1.8e-14 | 1.5e-14 | 1.5e-14 |
| `gfm_vctrl` | 300 | 6.2e-13 | 2.2e-13 | 3.7e-13 | 3.6e-13 |
| `gfl_pll1` | 500 | 3.5e-16 | 2.3e-13 | 2.8e-16 | 6.0e-13 |
| `gfl_delta` | 500 | 5.3e-14 | 2.0e-14 | 4.0e-14 | 1.9e-14 |
| `gfl_ip` | 500 | 0 | 0 | 4.2e-15 | 3.7e-14 |
| `gfl_iq` | 500 | 0 | 0 | 2.2e-15 | 1.3e-14 |
| `pmu_theta` | 1573 | 5.2e-14 | 2.5e-14 | 3.8e-14 | 2.3e-14 |
| `pmu_pi` | 1573 | 2.7e-16 | 1.3e-13 | 2.1e-16 | 3.2e-13 |
| `pmu_f` | 1573 | 2.1e-14 | 4.3e-16 | 2.1e-14 | 4.3e-16 |
| `cu_lat1` | 1573 | 6.6e-16 | 3.0e-13 | 6.5e-16 | 9.1e-13 |
| `cu_lat2` | 1573 | 5.6e-17 | 6.7e-16 | 4.2e-17 | 5.0e-16 |

The largest relative error is 8e-12 in the variant case and 4e-12 in the published case. Both
occur in `gfm_lat1`, whose values reach about 1e-6. The largest absolute error is 2e-12 in
`avr_lag`, whose values reach 3.6. `tools/compare_full.py` prints these tables from the full
logs.

## Solver

Simulink reports the solver ode3 with the step 0.01 s for the published setting
"FixedStepAuto". The port implements ode3 with the Butcher tableau and the order of sums of
Simulink. Two timings of the three Memory blocks were tried against the published case.

| Memory output in the minor time steps | Largest relative error after 1 s |
|---|---|
| the value of the last major time step (kept over the whole step) | 2e-12 |
| the new value, stored at the start of the step | 4e-2 |

The first timing agrees with Simulink. The port uses it.

## Tests

`tests/test_simulink.py` repeats the comparison with the smaller file
`tests/reference/wamces_reference.mat`. It holds all states at a few steps and three states of
each group over time. The published case keeps every tenth step of these traces. The variant
keeps every step. The tests require a relative error below 1e-10 for every group.

`tests/test_operating_point.py` checks without MATLAB that the undisturbed model starts at rest.
It also checks that the network gives the voltages of the power flow at the start.

## What the comparison does not cover

Simulink logs the continuous states only. The outputs of the Memory blocks and the bus voltages
were therefore compared through the states that depend on them. The disturbance acts from the
start in all cases, as in the published model. The model has no switching events or limits.

## Run times

Each case was simulated on one Windows machine that other jobs shared. Simulink took 152 s for the
published case in a first run and 421 s in a second run. The Python port took 92 s to 108 s.
These numbers show the order of magnitude only.
