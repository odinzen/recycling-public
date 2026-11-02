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


# House phase-field pastels (figure-color-where-strong): one color per field, black boundaries and labels.
FIELD_COLOR = {"beta": "#80b1d3", "two": "#ffffb3", "alpha": "#b3de69"}


def _figure(out) -> None:
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "legend.fontsize": 10})
    plt_, fig, axes = core.figure_axes(figsize=(6.5, 3.3), ncols=2)
    plt_.rcParams.update({"font.size": 10})
    for ax, el, tag in ((axes[0], "O", "(a)"), (axes[1], "N", "(b)")):
        rows = out[el]
        T = [r["T_K"] for r in rows]
        wa = [r["wt%_alpha"] for r in rows]
        wb = [r["wt%_beta"] for r in rows]
        xmax = 1.15 * max(wa)
        ax.fill_betweenx(T, 0, wb, color=FIELD_COLOR["beta"], linewidth=0)
        ax.fill_betweenx(T, wb, wa, color=FIELD_COLOR["two"], linewidth=0)
        ax.fill_betweenx(T, wa, xmax, color=FIELD_COLOR["alpha"], linewidth=0)
        ax.plot(wa, T, color="black", lw=1.4, marker="o", markerfacecolor="white", markeredgewidth=1.2, markersize=5, clip_on=False, zorder=3)
        ax.plot(wb, T, color="black", lw=1.4, ls=(0, (5, 2)), marker="s", markerfacecolor="black", markersize=5, clip_on=False, zorder=3)
        mid = len(T) // 2
        ax.text(0.5 * wb[-2], T[-2], "β", ha="center", va="center", fontsize=11)
        ax.text(0.5 * (wa[mid] + wb[mid]), T[mid], "α + β", ha="center", va="center", fontsize=10)
        ax.text(0.5 * (wa[1] + xmax), T[1] + 5, "α", ha="center", va="center", fontsize=11)
        ax.set_xlim(0, xmax)
        ax.set_ylim(min(T), max(T))
        ax.set_xlabel(f"{el} content (wt.%)")
        ax.text(0.04, 0.96, tag, transform=ax.transAxes, va="top", fontweight="bold")
        ax.tick_params(labelsize=10)
    axes[0].set_ylabel("Temperature (K)")
    fig.tight_layout()
    core.save(fig, "fig3_two_phase_field")
    plt_.close(fig)


if __name__ == "__main__":
    main()
