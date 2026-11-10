"""Input data of WAMCES: nine parameter tables and the power-flow solution.

The tables hold the values of the .mat files of the Zenodo record 10.5281/zenodo.18721887.
The script tools/make_data.py wrote them as CSV files. The power-flow solution is the
matrix ``bus_results`` of the file ``powerflow_sol.mat`` of that record. Every number is
stored with all the digits of its double value.
"""

import csv
from importlib import resources
from pathlib import Path

import numpy as np

# Short table names and the files that hold them.
TABLE_FILES = {
    "bus": "data_bus.csv",
    "lin": "data_lin.csv",
    "sm": "data_sm.csv",
    "avr": "data_avr.csv",
    "pss": "data_pss.csv",
    "tur": "data_tur.csv",
    "gfm": "data_gfm.csv",
    "gfl": "data_gfl.csv",
    "pmu": "data_pmu.csv",
}
POWERFLOW_FILE = "powerflow_sol.csv"


class Table(dict):
    """A table as a dict of column arrays. Columns can also be read as attributes."""

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name) from None

    def __setattr__(self, name, value):
        self[name] = value

    @property
    def nrows(self):
        return len(next(iter(self.values())))

    def copy(self):
        return Table({k: np.array(v, copy=True) for k, v in self.items()})


def _column(values):
    """Convert one CSV column to float, complex or str, in this order of preference."""
    try:
        return np.array([float(v) for v in values])
    except ValueError:
        pass
    try:
        # MATLAB writes complex numbers as "0-1.25i". Python expects "0-1.25j".
        return np.array([complex(v.replace("i", "j")) for v in values])
    except ValueError:
        return np.array(values, dtype=str)


def read_csv(path):
    """Read one CSV file with a header row into a Table."""
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    head, body = rows[0], rows[1:]
    return Table({name: _column([r[k] for r in body]) for k, name in enumerate(head)})


def data_dir():
    """The folder of the data files that ship with the package."""
    return Path(resources.files("wamces") / "data")


def load_tables(folder=None):
    """All nine tables and the power-flow solution as a dict of Tables.

    The key ``"pf"`` holds the power-flow solution with per-unit values on 100 MVA and
    angles in radians. A folder with the CSV files of the Zenodo record can replace the
    files of the package.
    """
    folder = Path(folder) if folder is not None else data_dir()
    tables = {key: read_csv(folder / name) for key, name in TABLE_FILES.items()}
    pf_path = folder / POWERFLOW_FILE
    if not pf_path.exists():
        pf_path = data_dir() / POWERFLOW_FILE
    tables["pf"] = read_csv(pf_path)
    return tables
