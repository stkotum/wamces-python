"""WAMCES in Python: the dynamic model of the continental European power system.

This package ports the Simulink model WAMCES of R. Musca (Zenodo, 2026,
doi:10.5281/zenodo.18721887) to Python. A short example::

    import wamces
    model = wamces.Model()              # 1000 MW load step at bus 2178
    result = model.simulate(t_end=20)   # fixed step 0.01 s, ode3
    result.f_coi                        # mean machine frequency in Hz
"""

from . import model
from .data import load_tables
from .dynamics import Dynamics, STATE_GROUPS
from .model import FN, SB, WN, Model
from .simulate import Result, simulate

__version__ = "1.0.0"
__all__ = ["Model", "simulate", "Result", "Dynamics", "load_tables", "STATE_GROUPS",
           "FN", "SB", "WN"]
