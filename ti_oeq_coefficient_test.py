"""Does the transus weight nitrogen at twice oxygen, as the oxygen equivalent O + 2N + 0.67C does?

The factor of two comes from hardness data. Here the alpha-stabilizing potency, the beta transus shift per wt% of
interstitial, is calculated from two independent assessments reduced to their metallic corner (Lindwall et al. 2018
for Ti-O, Ohtani and Hillert 1990 for Ti-N). Gate for each database: pure Ti at 1155 K.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import ti_transus_core as core

XS = [0.0, 0.002, 0.004, 0.006, 0.009, 0.012]
FIT_MAX_WT = 0.45


def curve(db, el):
    return [(core.wt_pct(x, el), core.transus(db, el, x)) for x in XS]


def slopes(rows):
    w, t = np.array(rows).T
    fit = float(np.polyfit(w[w <= FIT_MAX_WT], t[w <= FIT_MAX_WT], 1)[0])
    initial = float((t[1] - t[0]) / (w[1] - w[0]))
    return fit, initial


def main() -> None:
    out, curves = {}, {}
    for el in ("O", "N"):
        rows = curve(core.load(el), el)
        curves[el] = rows
        fit, initial = slopes(rows)
        out[el] = {"pure_Ti_K": round(rows[0][1], 2), "gate_pure_Ti_1155K": abs(rows[0][1] - 1155.0) < 0.5,
                   "slope_K_per_wt%": round(fit, 1), "initial_slope_K_per_wt%": round(initial, 1),
                   "rows": [{"wt%": round(w, 4), "transus_K": round(t, 2)} for w, t in rows]}
    ratio = out["N"]["slope_K_per_wt%"] / out["O"]["slope_K_per_wt%"]
    out["N_over_O_per_wt%"] = round(ratio, 2)
    out["N_over_O_per_wt%_dilute"] = round(out["N"]["initial_slope_K_per_wt%"] / out["O"]["initial_slope_K_per_wt%"], 2)
    out["N_over_O_per_atom"] = round(ratio * core.MASS["N"] / core.MASS["O"], 2)
    out["fit_range_wt%"] = f"0 to {FIT_MAX_WT}"
    out["O_eq_coefficient_for_N"] = 2.0
    Path(__file__).with_name("ti_oeq_coefficient_results.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    _figure(curves)


def _figure(curves) -> None:
    plt, fig, ax = core.figure_axes()
    for el, label, marker, fill, ls in (("O", "Ti-O", "o", "white", "-"), ("N", "Ti-N", "s", "black", (0, (5, 2)))):
        w, t = np.array(curves[el]).T
        ax.plot(w, t, color="black", lw=1.5, ls=ls, marker=marker, markerfacecolor=fill, markeredgewidth=1.3,
                markersize=7, label=label)
    ax.set_xlabel("Interstitial content (wt.%)")
    ax.set_ylabel("β transus (K)")
    ax.legend(loc="upper left", frameon=False)
    core.save(fig, "fig2_transus_O_N")
    plt.close(fig)


if __name__ == "__main__":
    main()
