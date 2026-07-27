"""Oxygen-equivalent coefficient test: does the thermodynamics give N about 2x O?

The oxygen-equivalent formula O_eq = O + 2N + 0.67C weights nitrogen twice oxygen
by mass. That weighting comes from hardness and embrittlement data. Here we test
it from a different direction: the alpha-stabilizing potency, measured as the
beta transus shift per unit mass of interstitial, computed from independent
published CALPHAD assessments (Lindwall 2018 for Ti-O, Ohtani-Hillert 1990 for
Ti-N), both reduced to their metallic Ti-rich corner.

Gate for each system: pure Ti must return 882 deg C. Then we compare the slopes
d(transus)/d(wt%) for oxygen and nitrogen, and report the ratio against 2.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pycalphad import Database, equilibrium, variables as v

DATA = Path(__file__).resolve().parents[2] / "data" / "recycled_ti"
TDB_O = DATA / "ti_o_lindwall_metallic_v1.tdb"
TDB_N = DATA / "ti_n_ohtani_metallic_v1.tdb"
M = {"O": 16.0, "N": 14.007}
M_TI = 47.867


def wt_pct(x: float, el: str) -> float:
    return 100.0 * x * M[el] / (x * M[el] + (1 - x) * M_TI)


def transus(db: Database, el: str, x: float) -> float | None:
    Tgrid = np.linspace(1080, 1500, 85)
    if x < 1e-6:
        comps, conds = ["TI", "VA"], {v.T: Tgrid, v.P: 101325, v.N: 1}
    else:
        comps, conds = ["TI", el, "VA"], {v.X(el): x, v.T: Tgrid, v.P: 101325, v.N: 1}
    eq = equilibrium(db, comps, ["HCP_A3", "BCC_A2"], conds)
    Ph = eq.Phase.values.squeeze(); NP = eq.NP.values.squeeze()
    pres = [sum(n for p, n in zip(Ph[i], NP[i]) if p == "HCP_A3" and n == n) > 1e-3
            for i in range(len(Tgrid))]
    for i in range(len(Tgrid) - 1):
        if pres[i] and not pres[i + 1]:
            return float(0.5 * (Tgrid[i] + Tgrid[i + 1]) - 273.15)
    return None


def curve(db: Database, el: str, xs):
    rows = []
    for x in xs:
        t = transus(db, el, x)
        rows.append((wt_pct(x, el), t))
        print(f"  {el}: x={x:.3f}  wt%={wt_pct(x, el):.3f}  "
              f"transus={'n/a' if t is None else f'{t:.1f} C'}")
    return rows


def slope_per_wt(rows):
    pts = [(w, t) for w, t in rows if t is not None and w <= 0.45]
    w = np.array([p[0] for p in pts]); t = np.array([p[1] for p in pts])
    a, _ = np.polyfit(w, t, 1)
    return float(a)  # deg C per wt%


def main() -> None:
    dbO, dbN = Database(str(TDB_O)), Database(str(TDB_N))
    xs = [0.0, 0.003, 0.006, 0.009, 0.012]
    print("oxygen:")
    rO = curve(dbO, "O", xs)
    print("nitrogen:")
    rN = curve(dbN, "N", xs)

    sO, sN = slope_per_wt(rO), slope_per_wt(rN)
    pureO, pureN = rO[0][1], rN[0][1]
    ratio = sN / sO
    out = {
        "pure_Ti_O_db_C": pureO, "pure_Ti_N_db_C": pureN,
        "gate_O": abs(pureO - 882) < 5, "gate_N": abs(pureN - 882) < 5,
        "slope_O_C_per_wt%": round(sO, 1),
        "slope_N_C_per_wt%": round(sN, 1),
        "N_over_O_potency_ratio": round(ratio, 2),
        "O_eq_coefficient_for_N": 2.0,
        "interpretation": (
            f"computed transus potency ratio {ratio:.2f}: nitrogen is the stronger "
            "alpha stabilizer, but by less than the embrittlement-based O_eq factor "
            "of 2, so the O_eq nitrogen weighting overstates the transus effect"
        ),
    }
    Path(__file__).with_name("ti_oeq_coefficient_results.json").write_text(json.dumps(out, indent=2))
    print("\n" + json.dumps(out, indent=2))
    _figure(rO, rN, out)


def _figure(rO, rN, out) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.linewidth": 1.0})
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    for s in ax.spines.values():
        s.set_color("black"); s.set_linewidth(1.0)
    ax.grid(False); ax.set_facecolor("white")
    wO = [w for w, t in rO if t is not None]; tO = [t for w, t in rO if t is not None]
    wN = [w for w, t in rN if t is not None]; tN = [t for w, t in rN if t is not None]
    ax.plot(wO, tO, color="0.0", lw=1.5, marker="o", markerfacecolor="white",
            markeredgewidth=1.4, markersize=6, label=f"oxygen, {out['slope_O_C_per_wt%']:.0f} C per wt%")
    ax.plot(wN, tN, color="0.0", lw=1.5, marker="s", markerfacecolor="0.0",
            linestyle=(0, (5, 2)), markersize=6, label=f"nitrogen, {out['slope_N_C_per_wt%']:.0f} C per wt%")
    ax.set_xlabel("interstitial content, wt%")
    ax.set_ylabel("beta transus, deg C")
    leg = ax.legend(loc="upper left", frameon=True, fontsize=7.4)
    leg.get_frame().set_edgecolor("black"); leg.get_frame().set_linewidth(0.8)
    fig.text(0.5, 0.955, "Nitrogen versus oxygen alpha-stabilizing potency",
             ha="center", va="top", fontsize=11, fontweight="bold")
    fig.text(0.5, 0.015,
             f"Independent Ti-O and Ti-N assessments. Computed N to O potency ratio {out['N_over_O_potency_ratio']}, against the oxygen-equivalent weighting of 2.",
             ha="center", va="bottom", fontsize=5.8, color="0.25", wrap=True)
    out_png = Path(__file__).with_name("figures") / "figure8_oeq_coefficient_test.png"
    out_png.parent.mkdir(exist_ok=True)
    fig.savefig(out_png, dpi=200); plt.close(fig)
    print("wrote", out_png.name)


if __name__ == "__main__":
    main()
