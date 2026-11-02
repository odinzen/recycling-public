"""Oxygen and nitrogen partitioning across the alpha + beta field, from the tie lines.

At each temperature in the two-phase field the coexisting alpha and beta compositions give the partition coefficient
k = (X/Ti in alpha) / (X/Ti in beta). The dilute-limit k from the same databases (ti_transus_core.dilute_partition)
is listed beside it: the two differ where the tie lines run to compositions at which the excess terms matter.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import ti_transus_core as core

TEMPS = np.arange(1160.0, 1241.0, 10.0)


def ratio(x):
    return x / (1 - x)


def main() -> None:
    out = {}
    for el in ("O", "N"):
        db = core.load(el)
        rows = []
        for T in TEMPS:
            tl = core.tie_line(db, el, float(T))
            if tl is None:
                continue
            xa, xb = tl
            rows.append({"T_K": float(T), "wt%_alpha": round(core.wt_pct(xa, el), 4),
                         "wt%_beta": round(core.wt_pct(xb, el), 4), "k_tie_line": round(ratio(xa) / ratio(xb), 3),
                         "k_dilute": round(core.dilute_partition(db, el, float(T)), 3)})
        out[el] = rows
        for r in rows:
            print(el, r)
    Path(__file__).with_name("ti_partition_results.json").write_text(json.dumps(out, indent=2))
    _figure(out)


def _figure(out) -> None:
    plt, fig, axes = core.figure_axes(figsize=(9.0, 4.0), ncols=2)
    for ax, el, name, tag in ((axes[0], "O", "Oxygen", "(a)"), (axes[1], "N", "Nitrogen", "(b)")):
        rows = out[el]
        T = [r["T_K"] for r in rows]
        wa = [r["wt%_alpha"] for r in rows]
        wb = [r["wt%_beta"] for r in rows]
        ax.fill_betweenx(T, wb, wa, hatch="////", facecolor="none", edgecolor="0.65", linewidth=0)
        ax.plot(wa, T, color="black", lw=1.5, marker="o", markerfacecolor="white", markeredgewidth=1.3,
                markersize=6, label="alpha boundary")
        ax.plot(wb, T, color="black", lw=1.5, ls=(0, (5, 2)), marker="s", markerfacecolor="black", markersize=6,
                label="beta boundary")
        ax.set_xlabel(f"{name} content (wt%)")
        ax.text(0.04, 0.95, tag, transform=ax.transAxes, va="top", fontweight="bold")
        ax.set_xlim(left=0)
    axes[0].set_ylabel("Temperature (K)")
    axes[1].legend(loc="lower right", frameon=False)
    fig.tight_layout()
    core.save(fig, "fig3_two_phase_field")
    plt.close(fig)


if __name__ == "__main__":
    main()
