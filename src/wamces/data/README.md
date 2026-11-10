# Data of WAMCES

The files in this folder come from the Zenodo record of WAMCES
([doi:10.5281/zenodo.18721887](https://doi.org/10.5281/zenodo.18721887), version 1, 2026) by
Rossano Musca. The record is licensed under CC BY 4.0 (see `LICENSE`). The model is described
by Musca, Ippolito and Riva Sanseverino in Electricity 7(2), 2026
([doi:10.3390/electricity7020028](https://doi.org/10.3390/electricity7020028)).

| File | Content | Source in the record |
|---|---|---|
| `data_bus.csv` | 3809 buses | `data_bus.mat` |
| `data_lin.csv` | 7343 lines | `data_lin.mat` |
| `data_sm.csv` | 618 synchronous machines | `data_sm.mat` |
| `data_avr.csv` | 618 voltage regulators | `data_avr.mat` |
| `data_pss.csv` | 618 power system stabilisers | `data_pss.mat` |
| `data_tur.csv` | 618 governors and turbines | `data_tur.mat` |
| `data_gfm.csv` | 300 grid-forming converters | `data_gfm.mat` |
| `data_gfl.csv` | 500 grid-following converters | `data_gfl.mat` |
| `data_pmu.csv` | 1573 phasor measurement units | `data_pmu.mat` |
| `powerflow_sol.csv` | Power-flow solution, per unit on 100 MVA, angles in rad | `powerflow_sol.mat` |

The files were written by `tools/make_data.py` from the .mat files of the record. They hold
the values of the .mat files exactly. The CSV files of the record (`raw_csv_data.zip`) are the
same tables with one difference. They round the column `Kf` of `data_gfm.csv` to 15 digits.
The file `powerflow_sol.csv` has no counterpart in `raw_csv_data.zip`.

The columns of `powerflow_sol.csv` are the columns of the matrix `bus_results`. The model uses
columns 3 to 8. The names of columns 1 and 2 (`ID`, `Type`) are given here. Column 2 equals
the column `Type` of `data_bus.csv`.
