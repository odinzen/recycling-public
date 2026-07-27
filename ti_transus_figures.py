"""Figures for the recycled Ti cases study (greyscale house style).

Renders the four manuscript figures. House style: pure greyscale, series
distinguished by marker, line texture, or diagonal hatching; full rectangular
box with four crisp black spines; white canvas with no gridlines; bold title at
the top of the canvas; a smaller source note at the bottom; axis labels
separating the variable name from its units with a comma.

Figure 2 uses literature calibration data (transus vs oxygen for Ti-6Al-4V).
Figures 1, 3, and 4 are schematics. No computed CALPHAD result is drawn; the
transus guide line in Figure 2 is the reported sensitivity, labeled as a guide,
not a database solve.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, FancyBboxPatch, Rectangle

OUTDIR = Path(__file__).with_name("figures")
OUTDIR.mkdir(exist_ok=True)

# Greyscale palette
INK = "0.0"
MIDGREY = "0.45"
LIGHTGREY = "0.82"
PANEL = "white"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.linewidth": 1.0,
        "axes.edgecolor": "black",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def _house_axes(ax) -> None:
    """Full box, four black spines, white canvas, no gridlines."""
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)
    ax.grid(False)
    ax.set_facecolor(PANEL)
    ax.tick_params(direction="out", color="black", length=4)


def _title_top(fig, text: str) -> None:
    fig.text(0.5, 0.955, text, ha="center", va="top", fontsize=11, fontweight="bold")


def _note_bottom(fig, text: str) -> None:
    fig.text(0.5, 0.015, text, ha="center", va="bottom", fontsize=6.3,
             color="0.25", wrap=True)


def _save(fig, name: str) -> Path:
    out = OUTDIR / name
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def figure1() -> Path:
    """PSPP qualification loop schematic."""
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.set_xticks([])
    ax.set_yticks([])
    _house_axes(ax)

    blocks = [
        (1.4, 3.4, "FEEDSTOCK\nchemistry, inclusions"),
        (5.0, 3.4, "GATES\nO_eq gate | inclusion gate"),
        (8.6, 3.4, "RESPONSE\ntransus surface,\nproperty ceiling"),
    ]
    for x, y, label in blocks:
        box = FancyBboxPatch(
            (x - 1.25, y - 0.7), 2.5, 1.4,
            boxstyle="round,pad=0.04", linewidth=1.2,
            edgecolor="black", facecolor=LIGHTGREY,
        )
        ax.add_patch(box)
        ax.text(x, y, label, ha="center", va="center", fontsize=8.2)

    for x0, x1 in [(2.65, 3.75), (6.25, 7.35)]:
        ax.add_patch(FancyArrow(x0, 3.4, x1 - x0, 0, width=0.02,
                                head_width=0.28, head_length=0.22,
                                length_includes_head=True, color=INK))
    # feedback arrow along the bottom
    ax.annotate("", xy=(1.4, 2.7), xytext=(8.6, 2.7),
                arrowprops=dict(arrowstyle="-|>", color=MIDGREY, lw=1.4,
                                connectionstyle="arc3,rad=-0.18"))
    ax.text(5.0, 1.55, "process iteration", ha="center", va="center",
            fontsize=8, color=MIDGREY, style="italic")

    _title_top(fig, "Qualification loop for recycled titanium")
    _note_bottom(fig, "Process, structure, property, performance loop after Olson 1997. Schematic.")
    fig.subplots_adjust(top=0.88, bottom=0.12, left=0.04, right=0.97)
    return _save(fig, "figure1_qualification_loop.png")


def figure2() -> Path:
    """Beta transus vs oxygen equivalent, Ti-6Al-4V; literature calibration."""
    # Literature calibration points (oxygen wt%, transus midpoint C, half-band C).
    # Nitrogen assumed negligible in these studies, so O_eq is taken as O.
    o_eq = [0.17, 0.20, 0.23]
    transus = [1002.5, 1007.5, 1017.5]
    band = [2.5, 2.5, 2.5]

    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    _house_axes(ax)

    # indicative slope guide: 5 to 10 C per 0.03 wt% O, drawn as the mid rate
    # (7.5 C / 0.03 = 250 C per wt%) through the central point. Labeled a guide.
    import numpy as np

    xg = np.linspace(0.14, 0.26, 50)
    yg = 1007.5 + 250.0 * (xg - 0.20)
    ax.plot(xg, yg, color=MIDGREY, linewidth=1.4, linestyle=(0, (6, 3)),
            label="indicative slope (guide), 7.5 C per 0.03 wt% O")

    ax.errorbar(o_eq, transus, yerr=band, fmt="s", color=INK, ecolor=INK,
                markersize=7, markerfacecolor="white", markeredgewidth=1.4,
                capsize=4, linewidth=0, label="literature calibration points")

    ax.set_xlabel("oxygen equivalent, wt%")
    ax.set_ylabel("beta transus, deg C")
    ax.set_xlim(0.14, 0.26)
    ax.set_ylim(995, 1025)
    leg = ax.legend(loc="upper left", frameon=True, fontsize=7.2)
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(0.8)

    _title_top(fig, "Beta transus against oxygen equivalent, Ti-6Al-4V")
    _note_bottom(fig, "Points: literature transus bands. Guide line is the reported sensitivity, not a database solve.")
    fig.subplots_adjust(top=0.88, bottom=0.16, left=0.15, right=0.96)
    return _save(fig, "figure2_transus_vs_oeq.png")


def figure3() -> Path:
    """Acceptance gate in the O_eq vs hydrogen plane."""
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    _house_axes(ax)

    # specification box: O_eq up to ~0.354, H up to 0.015
    oeq_spec, h_spec = 0.354, 0.015
    spec = Rectangle((0, 0), oeq_spec, h_spec, linewidth=1.4,
                     edgecolor="black", facecolor="white", zorder=2)
    ax.add_patch(spec)

    # embrittlement ceiling: O_eq above 0.25 is rejected; hatched region
    ceiling = Rectangle((0.25, 0), 0.13, h_spec, linewidth=0.0,
                        edgecolor="black", facecolor="none", hatch="////",
                        zorder=3)
    ax.add_patch(ceiling)
    ax.axvline(0.25, color=INK, linewidth=1.4, linestyle=(0, (4, 2)), zorder=4)

    # admissible region annotation
    ax.text(0.125, 0.0075, "admissible", ha="center", va="center", fontsize=9.5)
    ax.text(0.315, 0.0075, "embrittled\n(rejected)", ha="center", va="center",
            fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor="none"))
    # ceiling label placed in the white admissible region, left of the line
    ax.text(0.243, 0.0165, "O_eq ceiling, 0.25 wt%", ha="right", va="center",
            fontsize=7.4, color=INK)

    ax.set_xlabel("oxygen equivalent, wt%")
    ax.set_ylabel("hydrogen, wt%")
    ax.set_xlim(0, 0.40)
    ax.set_ylim(0, 0.021)

    _title_top(fig, "Acceptance gate, oxygen equivalent against hydrogen")
    _note_bottom(fig, "Box: grade 5 maxima. Hatched: composite O_eq ceiling cutting inside the box.")
    fig.subplots_adjust(top=0.88, bottom=0.16, left=0.16, right=0.96)
    return _save(fig, "figure3_acceptance_gate.png")


def figure4() -> Path:
    """Cold hearth inclusion removal schematic."""
    import numpy as np

    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.set_xticks([])
    ax.set_yticks([])
    _house_axes(ax)

    # melt pool
    pool = FancyBboxPatch((1.0, 1.0), 8.0, 2.6, boxstyle="round,pad=0.02",
                          linewidth=1.2, edgecolor="black", facecolor=LIGHTGREY)
    ax.add_patch(pool)
    ax.text(5.0, 3.35, "cold hearth melt pool", ha="center", va="center",
            fontsize=8.5, style="italic")

    # flow direction
    ax.add_patch(FancyArrow(1.3, 3.0, 7.3, 0, width=0.015, head_width=0.22,
                            head_length=0.2, length_includes_head=True,
                            color=MIDGREY))
    ax.text(1.2, 4.2, "feed", ha="left", fontsize=8)
    ax.text(9.0, 4.2, "ingot", ha="right", fontsize=8)

    # low density inclusions dissolving (shrinking open circles)
    rng = np.linspace(2.2, 7.5, 5)
    for i, x in enumerate(rng):
        r = 0.22 * (1 - i / 6.0)
        ax.add_patch(plt.Circle((x, 2.6), r, fill=False, edgecolor=INK, lw=1.2))
    ax.text(4.8, 2.05, "low density inclusions dissolve", ha="center",
            fontsize=7.4)

    # high density inclusion settling (filled square with down arrow)
    ax.add_patch(Rectangle((6.7, 2.4), 0.26, 0.26, facecolor=INK))
    ax.add_patch(FancyArrow(6.83, 2.35, 0, -0.9, width=0.03, head_width=0.18,
                            head_length=0.18, length_includes_head=True,
                            color=INK))
    ax.text(7.15, 1.55, "high density\ninclusion settles", ha="left",
            va="center", fontsize=7.4)

    _title_top(fig, "Inclusion removal in cold hearth melting")
    _note_bottom(fig, "Residence time dissolves low density inclusions; density settles high density ones. Schematic.")
    fig.subplots_adjust(top=0.88, bottom=0.1, left=0.04, right=0.97)
    return _save(fig, "figure4_inclusion_removal.png")


def main() -> list[Path]:
    paths = [figure1(), figure2(), figure3(), figure4()]
    for p in paths:
        print("wrote", p.name)
    return paths


if __name__ == "__main__":
    main()
