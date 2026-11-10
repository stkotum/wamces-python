"""The effect of the wide-area gain Kw of the grid-forming converters.

Each grid-forming converter receives the frequency that the central unit measures from
all PMUs, after the latency Td. The gain Kw weights the difference between that frequency
and the converter's own speed. The published script sets Kw to zero. This example runs
the published load step for several values of Kw. It prints three numbers per run and
plots the mean frequency and the spread of the machine frequencies.

Run: python examples/wide_area_gain.py [figure.png]
The four runs take five to ten minutes.
"""

import sys

import matplotlib.pyplot as plt
import numpy as np

import wamces

GAINS = (0.0, 10.0, 20.0, 40.0)
T_END = 20.0
# An ordinal blue ramp: a darker line is a larger gain.
COLORS = ("#6da7ec", "#2a78d6", "#1c5cab", "#0d366b")
SURFACE, INK, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"


def run(kw):
    result = wamces.Model(kw=kw).simulate(t_end=T_END, save_every=1, save=["sm_domega"])
    f = result.f_sm
    dev = f - wamces.FN * (1.0 + result.dfsys)[:, None]  # each machine against the mean
    return result, {
        "nadir": result.f_coi.min(),
        "spread_max": np.ptp(f, axis=1).max(),
        "dev_rms": np.sqrt(np.mean(dev**2)),
    }


def main(path=None):
    runs = [run(kw) for kw in GAINS]
    print("| Kw | lowest mean frequency (Hz) | largest spread (mHz) | RMS deviation from the mean (mHz) |")
    print("|---|---|---|---|")
    for kw, (_, m) in zip(GAINS, runs):
        print("| %g | %.4f | %.2f | %.3f |" % (kw, m["nadir"], 1e3 * m["spread_max"], 1e3 * m["dev_rms"]))

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
                         "grid.color": GRID, "grid.linewidth": 0.8})
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(8, 6), facecolor=SURFACE)
    for (result, _), kw, c in zip(runs, GAINS, COLORS):
        ax[0].plot(result.t_signals, result.f_coi, color=c, lw=2, label="Kw = %g" % kw)
        ax[1].plot(result.t, 1e3 * np.ptp(result.f_sm, axis=1), color=c, lw=2)
    ax[0].set_ylabel("mean machine frequency (Hz)")
    ax[1].set_ylabel("spread of the machine\nfrequencies (mHz)")
    ax[1].set_xlabel("time (s)")
    ax[0].legend(frameon=False, labelcolor=INK)
    for a in ax:
        a.set_facecolor(SURFACE)
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
    fig.align_ylabels(ax)
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150, facecolor=SURFACE)
    else:
        plt.show()


if __name__ == "__main__":
    main(*sys.argv[1:2])
