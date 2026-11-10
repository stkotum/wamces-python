"""The published case: a load step of 1000 MW at bus 2178, simulated for 50 s.

Run: python examples/published_case.py [figure.png]
The top panel shows the mean frequency of the synchronous machines and the frequency
that the central unit measures. The bottom panel shows the frequency of every machine.
The run takes one to three minutes.
"""

import sys

import matplotlib.pyplot as plt

import wamces

SURFACE, INK, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, MACHINES = "#2a78d6", "#eb6834", "#9ec5f4"


def main(path=None):
    result = wamces.Model().simulate(t_end=50.0, save_every=5, save=["sm_domega"])

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
                         "grid.color": GRID, "grid.linewidth": 0.8})
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(8, 6), facecolor=SURFACE)
    ax[0].plot(result.t_signals, result.f_coi, color=BLUE, lw=2, label="mean of the machines")
    ax[0].plot(result.t_signals, result.f_wide_area, color=ORANGE, lw=2, ls="--",
               label="central unit (PMUs)")
    ax[0].set_ylabel("frequency (Hz)")
    ax[0].legend(frameon=False, labelcolor=INK)
    ax[1].plot(result.t, result.f_sm, color=MACHINES, lw=0.6)
    ax[1].plot(result.t_signals, result.f_coi, color=BLUE, lw=2)
    ax[1].set_ylabel("frequency of the\n618 machines (Hz)")
    ax[1].set_xlabel("time (s)")
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
