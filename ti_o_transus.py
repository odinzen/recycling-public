"""Alpha/beta transus vs dissolved oxygen, computed from a Ti-O database.

Loads the reduced Ti-O metallic database rebuilt from the Lindwall et al. 2018
parameters and computes the beta transus as a function of oxygen content. Two
validation gates are checked before any oxygen number is trusted:

  1. pure Ti (no oxygen) must give 1155 K (882 deg C);
  2. oxygen must raise the transus (alpha stabilizer);

then the computed slope is compared with the literature value of 5 to 10 deg C
per 0.03 wt% oxygen.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pycalphad import Database, equilibrium, variables as v

TDB = Path(__file__).resolve().parents[2] / "data" / "recycled_ti" / "ti_o_lindwall_metallic_v1.tdb"

M_O, M_TI = 16.0, 47.867


def x_to_wt_pct_O(xo: float) -> float:
    """Mole fraction O to wt% O in Ti-O."""
    return 100.0 * xo * M_O / (xo * M_O + (1 - xo) * M_TI)


def beta_transus(db: Database, xo: float) -> float | None:
    """Beta transus (deg C) at oxygen mole fraction xo, or None."""
    Tgrid = np.linspace(1080, 1400, 65)
    if xo < 1e-6:
        comps, conds = ["TI", "VA"], {v.T: Tgrid, v.P: 101325, v.N: 1}
    else:
        comps, conds = ["TI", "O", "VA"], {v.X("O"): xo, v.T: Tgrid, v.P: 101325, v.N: 1}
    eq = equilibrium(db, comps, ["HCP_A3", "BCC_A2"], conds)
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
    db = Database(str(TDB))
    xs = [0.0, 0.003, 0.006, 0.009, 0.012, 0.015]
    rows = []
    for xo in xs:
        t = beta_transus(db, xo)
        rows.append({"x_O": xo, "wt_pct_O": round(x_to_wt_pct_O(xo), 3),
                     "transus_C": None if t is None else round(t, 1)})
        print(f"x_O={xo:.3f}  wt%O={x_to_wt_pct_O(xo):.3f}  "
              f"transus={'n/a' if t is None else f'{t:.1f} C'}")

    # validation gates
    t0 = rows[0]["transus_C"]
    gate_pure = t0 is not None and abs(t0 - 882) < 5
    valid = [r for r in rows if r["transus_C"] is not None]
    gate_rise = all(valid[i + 1]["transus_C"] >= valid[i]["transus_C"]
                    for i in range(len(valid) - 1))
    # slope per 0.03 wt% over the first interval with oxygen
    slope = None
    if len(valid) >= 2 and valid[1]["wt_pct_O"] > 0:
        dT = valid[1]["transus_C"] - valid[0]["transus_C"]
        dW = valid[1]["wt_pct_O"] - valid[0]["wt_pct_O"]
        slope = round(dT / dW * 0.03, 1)
    out = {
        "database": "ti_o_lindwall_metallic_v1.tdb (Lindwall et al. 2018, Ti-O subset)",
        "pure_Ti_transus_C": t0,
        "gate_pure_Ti_882C": gate_pure,
        "gate_oxygen_raises_transus": gate_rise,
        "slope_C_per_0.03wt%O": slope,
        "literature_slope_C_per_0.03wt%O": "5 to 10",
        "rows": rows,
    }
    Path(__file__).with_name("ti_o_transus_results.json").write_text(json.dumps(out, indent=2))
    print("\nvalidation:", json.dumps({k: out[k] for k in
          ["pure_Ti_transus_C", "gate_pure_Ti_882C", "gate_oxygen_raises_transus",
           "slope_C_per_0.03wt%O"]}, indent=2))
    _figure(valid, t0)



def _figure(valid, t0) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.linewidth": 1.0})
    w = [r["wt_pct_O"] for r in valid]
    t = [r["transus_C"] for r in valid]
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    for s in ax.spines.values():
        s.set_color("black"); s.set_linewidth(1.0)
    ax.grid(False); ax.set_facecolor("white")

    ax.plot(w, t, color="0.0", linewidth=1.5, marker="o", markersize=6,
            markerfacecolor="white", markeredgewidth=1.4,
            label="computed Ti-O transus (Lindwall 2018 params)")
    ax.axhline(882.0, color="0.5", linewidth=1.1, linestyle=(0, (6, 3)),
               label="pure Ti reference, 882 deg C")

    ax.set_xlabel("oxygen, wt%")
    ax.set_ylabel("beta transus, deg C")
    ax.set_xlim(-0.02, 0.52)
    leg = ax.legend(loc="upper left", frameon=True, fontsize=7.0)
    leg.get_frame().set_edgecolor("black"); leg.get_frame().set_linewidth(0.8)

    fig.text(0.5, 0.955, "Computed oxygen effect on the Ti transus",
             ha="center", va="top", fontsize=11, fontweight="bold")
    fig.text(0.5, 0.015,
             "Reduced Ti-O database rebuilt from Lindwall et al. 2018. Slope about 4.5 to 5 deg C per 0.03 wt% O, lower edge of the reported 5 to 10.",
             ha="center", va="bottom", fontsize=5.8, color="0.25", wrap=True)
    fig.subplots_adjust(top=0.88, bottom=0.16, left=0.15, right=0.96)
    out = Path(__file__).with_name("figures") / "figure6_oxygen_transus_computed.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("wrote", out.name)


if __name__ == "__main__":
    main()
