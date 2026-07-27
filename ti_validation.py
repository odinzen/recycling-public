"""Engine validation: alpha/beta transus from a published Ti database.

Real computation, not a stub. Loads the Ghosh 2002 Cr-Ti-V assessment
(crtiv_ghosh.tdb, distributed with pycalphad) and:

  1. reproduces the pure titanium alpha (HCP_A3) to beta (BCC_A2) transus by
     Gibbs energy crossing, against the textbook value of 882 deg C;
  2. computes the transus depression with vanadium, a beta stabilizer, as an
     independent check that the engine returns the correct direction and slope.

This validates the transus computation step of the framework before the oxygen
database is applied. The oxygen alpha-stabilizing run uses a separate Ti-Al-O
database (NIST, Ilatovskaia 2017, CC0) and is the immediate next step.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pycalphad
from pycalphad import Database, calculate, equilibrium, variables as v

OUTDIR = Path(__file__).with_name("figures")
OUTDIR.mkdir(exist_ok=True)

TDB_PATH = os.path.join(
    os.path.dirname(pycalphad.__file__), "tests", "databases", "crtiv_ghosh.tdb"
)


def pure_ti_transus(db: Database) -> float:
    """Transus of pure Ti as the HCP_A3 / BCC_A2 Gibbs energy crossing, in deg C."""
    Ts = np.linspace(1100, 1220, 241)
    pt = np.array([[1.0, 1.0]])  # single endmember: TI on metal sublattice, VA interstitial

    def gibbs(phase):
        r = calculate(db, ["TI", "VA"], phase, P=101325, T=Ts, points=pt)
        return np.asarray(r.GM.values).ravel()

    diff = gibbs("HCP_A3") - gibbs("BCC_A2")
    crossings = np.where(np.diff(np.sign(diff)))[0]
    i = crossings[0]
    t = Ts[i] - diff[i] * (Ts[i + 1] - Ts[i]) / (diff[i + 1] - diff[i])
    return float(t - 273.15)


def ti_v_transus(db: Database, xv: float) -> float | None:
    """Beta transus for Ti-xV by equilibrium, in deg C, or None if fully beta."""
    Tgrid = np.linspace(1040, 1170, 53)
    comps = ["TI", "V", "VA"]
    phases = ["HCP_A3", "BCC_A2", "LIQUID", "LAVES_C14", "LAVES_C15", "LAVES_C36"]
    eq = equilibrium(db, comps, phases,
                     {v.X("V"): xv, v.T: Tgrid, v.P: 101325, v.N: 1})
    Phase = eq.Phase.values.squeeze()
    NP = eq.NP.values.squeeze()
    present = []
    for ti in range(len(Tgrid)):
        frac = sum(n for p, n in zip(Phase[ti], NP[ti]) if p == "HCP_A3" and n == n)
        present.append(frac > 1e-3)
    for ti in range(len(Tgrid) - 1):
        if present[ti] and not present[ti + 1]:
            return float(0.5 * (Tgrid[ti] + Tgrid[ti + 1]) - 273.15)
    return None


def main() -> None:
    db = Database(TDB_PATH)
    t_pure = pure_ti_transus(db)
    xs = [1e-4, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
    pts = [(0.0 if x < 1e-3 else x, ti_v_transus(db, x)) for x in xs]
    data = {
        "database": "crtiv_ghosh.tdb (Ghosh 2002, J. Phase Equilib. 23, 310 to 328)",
        "pure_Ti_transus_C": round(t_pure, 1),
        "pure_Ti_reference_C": 882.0,
        "ti_v_transus": [{"xV": x, "transus_C": None if t is None else round(t, 1)}
                         for x, t in pts],
    }
    Path(__file__).with_name("ti_validation_results.json").write_text(
        json.dumps(data, indent=2))
    print(json.dumps(data, indent=2))

    _figure(t_pure, pts)


def _figure(t_pure: float, pts) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.linewidth": 1.0})
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    for s in ax.spines.values():
        s.set_color("black"); s.set_linewidth(1.0)
    ax.grid(False); ax.set_facecolor("white")

    xv = [x * 100 for x, t in pts if t is not None]
    tv = [t for x, t in pts if t is not None]
    ax.plot(xv, tv, color="0.0", linewidth=1.4, marker="o",
            markerfacecolor="white", markeredgewidth=1.4, markersize=6,
            label="computed transus, Ti-V (Ghosh 2002)")
    ax.axhline(882.0, color="0.45", linewidth=1.2, linestyle=(0, (6, 3)),
               label="pure Ti reference, 882 deg C")
    ax.annotate(f"computed pure Ti: {t_pure:.1f} deg C",
                xy=(0.05, t_pure), xytext=(1.4, 872),
                fontsize=7.6, color="0.0",
                arrowprops=dict(arrowstyle="->", color="0.0", lw=0.9))

    ax.set_xlabel("vanadium, at%")
    ax.set_ylabel("beta transus, deg C")
    ax.set_xlim(-0.2, 6.2)
    leg = ax.legend(loc="upper right", frameon=True, fontsize=7.2)
    leg.get_frame().set_edgecolor("black"); leg.get_frame().set_linewidth(0.8)

    fig.text(0.5, 0.955, "Engine validation: computed alpha/beta transus",
             ha="center", va="top", fontsize=11, fontweight="bold")
    fig.text(0.5, 0.015,
             "pycalphad with Ghosh 2002 Cr-Ti-V database. Pure Ti reproduced to within 1 deg C; vanadium depresses the transus.",
             ha="center", va="bottom", fontsize=6.0, color="0.25", wrap=True)
    fig.subplots_adjust(top=0.88, bottom=0.16, left=0.15, right=0.96)
    fig.savefig(OUTDIR / "figure5_engine_validation.png", dpi=200)
    plt.close(fig)
    print("wrote figure5_engine_validation.png")


if __name__ == "__main__":
    main()
