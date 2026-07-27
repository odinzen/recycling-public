"""Oxygen partitioning and the alpha+beta two-phase field in Ti-O.

A phase-equilibria result the single transus line does not give: at each
temperature in the two-phase field, equilibrium fixes the oxygen content of the
coexisting alpha (HCP_A3) and beta (BCC_A2) phases. Their ratio is the oxygen
partition coefficient k_O = wt%O(alpha) / wt%O(beta), which quantifies how
strongly oxygen concentrates in alpha during heat treatment of recycled Ti.

Computed from data/recycled_ti/ti_o_lindwall_metallic_v1.tdb.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pycalphad import Database, equilibrium, variables as v

TDB = Path(__file__).resolve().parents[2] / "data" / "recycled_ti" / "ti_o_lindwall_metallic_v1.tdb"
M_O, M_TI = 16.0, 47.867


def wt_pct_O(xo: float) -> float:
    return 100.0 * xo * M_O / (xo * M_O + (1 - xo) * M_TI)


def tie_line(db: Database, T: float) -> dict | None:
    """At temperature T (K), return oxygen content of coexisting alpha and beta."""
    # scan overall oxygen; in the two-phase field both phases appear with
    # constant compositions (the tie-line endpoints).
    xo_grid = np.linspace(0.0005, 0.020, 40)
    eq = equilibrium(db, ["TI", "O", "VA"], ["HCP_A3", "BCC_A2"],
                     {v.X("O"): xo_grid, v.T: T, v.P: 101325, v.N: 1})
    Phase = eq.Phase.values.squeeze()
    XO = eq.X.sel(component="O").values.squeeze()  # (xo, vertex)
    NP = eq.NP.values.squeeze()
    for i in range(len(xo_grid)):
        ph = list(Phase[i]); npf = list(NP[i]); xo = list(XO[i])
        present = {p: (n, x) for p, n, x in zip(ph, npf, xo) if isinstance(p, str) and p}
        if "HCP_A3" in present and "BCC_A2" in present \
                and present["HCP_A3"][0] > 1e-3 and present["BCC_A2"][0] > 1e-3:
            xa = float(present["HCP_A3"][1]); xb = float(present["BCC_A2"][1])
            wa, wb = wt_pct_O(xa), wt_pct_O(xb)
            return {"T_C": round(T - 273.15, 1),
                    "wtO_alpha": round(wa, 4), "wtO_beta": round(wb, 4),
                    "k_O": round(wa / wb, 2) if wb > 0 else None}
    return None


def main() -> None:
    db = Database(str(TDB))
    rows = []
    for T in np.arange(1133, 1216, 10.0):  # ~860 to 940 C
        tl = tie_line(db, float(T))
        if tl:
            rows.append(tl)
            print(f"T={tl['T_C']} C  alpha={tl['wtO_alpha']} wt%O  "
                  f"beta={tl['wtO_beta']} wt%O  k_O={tl['k_O']}")
    out = {"database": "ti_o_lindwall_metallic_v1.tdb (Lindwall 2018 Ti-O subset)",
           "note": "k_O = wt%O(alpha)/wt%O(beta); oxygen concentrates in alpha (k_O>1).",
           "rows": rows}
    Path(__file__).with_name("ti_o_partition_results.json").write_text(json.dumps(out, indent=2))
    if rows:
        ks = [r["k_O"] for r in rows if r["k_O"]]
        print(f"\nk_O range: {min(ks)} to {max(ks)}  (oxygen prefers alpha)")
    _figure(rows)


def _figure(rows) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not rows:
        return
    Tc = [r["T_C"] for r in rows]
    wa = [r["wtO_alpha"] for r in rows]
    wb = [r["wtO_beta"] for r in rows]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.linewidth": 1.0})
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    for s in ax.spines.values():
        s.set_color("black"); s.set_linewidth(1.0)
    ax.grid(False); ax.set_facecolor("white")
    ax.plot(wa, Tc, color="0.0", linewidth=1.5, marker="o", markersize=5,
            markerfacecolor="white", markeredgewidth=1.3, label="alpha boundary")
    ax.plot(wb, Tc, color="0.0", linewidth=1.5, marker="s", markersize=5,
            markerfacecolor="0.0", linestyle=(0, (5, 2)), label="beta boundary")
    ax.fill_betweenx(Tc, wb, wa, hatch="////", facecolor="none", edgecolor="0.6", linewidth=0)
    ax.set_xlabel("oxygen, wt%")
    ax.set_ylabel("temperature, deg C")
    leg = ax.legend(loc="lower right", frameon=True, fontsize=7.4)
    leg.get_frame().set_edgecolor("black"); leg.get_frame().set_linewidth(0.8)
    fig.text(0.5, 0.955, "Oxygen partitioning across the alpha + beta field",
             ha="center", va="top", fontsize=11, fontweight="bold")
    fig.text(0.5, 0.015,
             "Computed tie-lines, Ti-O. Oxygen concentrates in alpha; the two phase field between the curves widens with oxygen.",
             ha="center", va="bottom", fontsize=6.0, color="0.25", wrap=True)
    fig.subplots_adjust(top=0.88, bottom=0.16, left=0.15, right=0.96)
    out = Path(__file__).with_name("figures") / "figure7_oxygen_partition.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("wrote", out.name)


if __name__ == "__main__":
    main()
