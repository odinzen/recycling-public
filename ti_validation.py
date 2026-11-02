"""Validation of the transus calculation on a published Ti database (Ghosh 2002 Cr-Ti-V, distributed with pycalphad).

Two checks before any interstitial result is used: pure Ti must transform at 1155 K, found here as the HCP_A3 /
BCC_A2 Gibbs energy crossing, and vanadium, a beta stabilizer, must lower the transus, found by bisection on the
temperature at which the last alpha disappears.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pycalphad
from pycalphad import Database, calculate, equilibrium, variables as v

import ti_transus_core as core

TDB_PATH = os.path.join(os.path.dirname(pycalphad.__file__), "tests", "databases", "crtiv_ghosh.tdb")
PHASES = ["HCP_A3", "BCC_A2", "LIQUID", "LAVES_C14", "LAVES_C15", "LAVES_C36"]
M_V = 50.942
WT_V = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]


def pure_ti_transus(db: Database) -> float:
    Ts = np.linspace(1100, 1220, 241)
    pt = np.array([[1.0, 1.0]])
    diff = [np.asarray(calculate(db, ["TI", "VA"], ph, P=101325, T=Ts, points=pt).GM.values).ravel()
            for ph in ("HCP_A3", "BCC_A2")]
    diff = diff[0] - diff[1]
    i = np.where(np.diff(np.sign(diff)))[0][0]
    return float(Ts[i] - diff[i] * (Ts[i + 1] - Ts[i]) / (diff[i + 1] - diff[i]))


def x_from_wt(wt: float) -> float:
    n_v, n_ti = wt / M_V, (100 - wt) / core.MASS["TI"]
    return n_v / (n_v + n_ti)


def hcp_fraction(db, x, T):
    eq = equilibrium(db, ["TI", "V", "VA"], PHASES, {v.X("V"): x, v.T: T, v.P: 101325, v.N: 1})
    ph, n = eq.Phase.values.squeeze(), eq.NP.values.squeeze()
    return float(sum(nn for p, nn in zip(ph, n) if p == "HCP_A3" and nn == nn))


def ti_v_transus(db, wt, lo=950.0, hi=1200.0, tol=0.01) -> float:
    x = max(x_from_wt(wt), 1e-6)
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if hcp_fraction(db, x, mid) > 1e-6 else (lo, mid)
    return 0.5 * (lo + hi)


def main() -> None:
    db = Database(TDB_PATH)
    t_pure = pure_ti_transus(db)
    rows = [{"wt%_V": w, "transus_K": round(ti_v_transus(db, w), 2)} for w in WT_V]
    data = {"database": "crtiv_ghosh.tdb (Ghosh 2002, J. Phase Equilib. 23, 310 to 328)",
            "pure_Ti_crossing_K": round(t_pure, 2), "pure_Ti_reference_K": 1155.0, "ti_v_transus": rows}
    Path(__file__).with_name("ti_validation_results.json").write_text(json.dumps(data, indent=2))
    print(json.dumps(data, indent=2))
    plt, fig, ax = core.figure_axes()
    ax.plot([r["wt%_V"] for r in rows], [r["transus_K"] for r in rows], color="black", lw=1.5, marker="o",
            markerfacecolor="white", markeredgewidth=1.3, markersize=7, label="Ti-V, calculated")
    ax.axhline(1155.0, color="0.45", lw=1.2, ls=(0, (6, 3)), label="pure Ti, 1155 K")
    ax.set_xlabel("V content (wt.%)")
    ax.set_ylabel("β transus (K)")
    ax.legend(loc="lower left", frameon=False)
    core.save(fig, "fig1_validation_ti_v")
    plt.close(fig)


if __name__ == "__main__":
    main()
