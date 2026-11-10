"""Write the CSV files of the package from the published .mat files.

Usage:
    python tools/make_data.py RAW_TABLES_MAT POWERFLOW_SOL_MAT OUTPUT_FOLDER

``RAW_TABLES_MAT`` is the output of ``tools/export_tables.m``. ``POWERFLOW_SOL_MAT`` is
the file ``powerflow_sol.mat`` of the Zenodo record. Every number is written with the
shortest text that reads back to the same double. The files then hold the values of the
.mat files exactly. The published CSV files round some values to 15 digits.
"""

import csv
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio

TABLES = {"bus": "Tbus", "lin": "Tlin", "sm": "Tsm", "avr": "Tavr", "pss": "Tpss",
          "tur": "Ttur", "gfm": "Tgfm", "gfl": "Tgfl", "pmu": "Tpmu"}
POWERFLOW_COLUMNS = ["ID", "Type", "Vmag", "Vang", "Pgen", "Qgen", "Ploa", "Qloa"]


def fmt(v):
    """Shortest text of a value. Complex numbers take the MATLAB form 0-1.25i."""
    if isinstance(v, str):
        return v
    if isinstance(v, (complex, np.complexfloating)):
        sign = "-" if np.signbit(v.imag) else "+"
        return "%s%s%si" % (fmt(v.real), sign, fmt(abs(v.imag)))
    f = float(v)
    if f == int(f) and abs(f) < 1e15:
        return str(int(f))
    return repr(f)


def write(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        for row in rows:
            w.writerow([fmt(x) for x in row])


def main(raw_path, pf_path, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    raw = sio.loadmat(raw_path, squeeze_me=True, struct_as_record=False)
    for key, var in TABLES.items():
        cols = list(raw[var + "_columns"])
        data = [np.atleast_1d(getattr(raw[var], c)) for c in cols]
        write(out / ("data_%s.csv" % key), cols, zip(*data))
    pf = sio.loadmat(pf_path)["bus_results"]
    write(out / "powerflow_sol.csv", POWERFLOW_COLUMNS, pf)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
