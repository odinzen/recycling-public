"""Enthalpy and entropy of moving O or N from beta to alpha Ti, and how they set the transus potency.

In the dilute limit the two reduced databases give, per mole of interstitial,

    ln k = ln(a_alpha / a_beta) + dS_nc / R - dH / (R T),

with a the interstitial sites per Ti (0.5 in alpha, 3 in beta), so the configurational term is the same for O and N
and cancels from their comparison. The transus shift follows from the same k (equal Ti chemical potentials,
dilute solution):

    dT / r_beta = R T0^2 (k - 1) / dH_Ti,

r_beta the X/Ti ratio of the beta at its transus and dH_Ti the alpha-to-beta enthalpy of pure Ti. The script checks
that relation against the calculated transus and evaluates the N:O potency with the two entropies made equal.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import ti_transus_core as core

T_REF = 1180.0  # middle of the calculated two-phase field
T0 = 1155.0


def per_wt_slope(k, el, dH_ti):
    """Dilute transus slope, K per wt%, from k."""
    per_ratio = core.R * T0**2 * (k - 1) / dH_ti  # K per unit X/Ti
    return per_ratio * core.MASS["TI"] / (100 * core.MASS[el])


def main() -> None:
    dbs = {el: core.load(el) for el in ("O", "N")}
    dH_ti = core.ti_transformation_enthalpy(dbs["O"], T0)
    split = {el: core.transfer_split(dbs[el], el, T_REF) for el in dbs}
    k0 = {el: core.dilute_partition(dbs[el], el, T0) for el in dbs}
    slope = {el: per_wt_slope(k0[el], el, dH_ti) for el in dbs}
    out = {"T_ref_K": T_REF, "dH_Ti_alpha_to_beta_J_mol": round(dH_ti, 1), "split": {}, "dilute_k_at_1155K": {}}
    for el in dbs:
        s = split[el]
        out["split"][el] = {key: round(float(val), 3) for key, val in s.items()}
        out["dilute_k_at_1155K"][el] = round(k0[el], 3)
    out["dilute_slope_from_k_K_per_wt%"] = {el: round(slope[el], 1) for el in dbs}
    out["N_over_O_per_wt%_from_k"] = round(slope["N"] / slope["O"], 2)

    # The same comparison with the two interstitials given one non-configurational entropy, each choice in turn.
    def k_with(el, s_nc):
        s = split[el]
        return float(np.exp(s["dS_config"] / core.R + s_nc / core.R - s["dH"] / (core.R * T0)))
    for basis in ("O", "N"):
        s_nc = split[basis]["dS_nonconfig"]
        ratio = per_wt_slope(k_with("N", s_nc), "N", dH_ti) / per_wt_slope(k_with("O", s_nc), "O", dH_ti)
        out[f"N_over_O_per_wt%_if_both_had_{basis}_entropy"] = round(ratio, 2)
    Path(__file__).with_name("ti_transfer_entropy_results.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    _figure(dbs)


def _figure(dbs) -> None:
    part = json.loads(Path(__file__).with_name("ti_partition_results.json").read_text())
    plt, fig, ax = core.figure_axes()
    T = np.linspace(1150, 1250, 41)
    for el, name, marker, fill, ls in (("O", "O", "o", "white", "-"), ("N", "N", "s", "black", (0, (5, 2)))):
        ax.plot(1000 / T, [np.log(core.dilute_partition(dbs[el], el, t)) for t in T], color="black", lw=1.5, ls=ls,
                label=f"{name}, dilute limit")
        rows = part[el]
        ax.plot([1000 / r["T_K"] for r in rows], [np.log(r["k_tie_line"]) for r in rows], ls="none", marker=marker,
                markerfacecolor=fill, markeredgecolor="black", markeredgewidth=1.3, markersize=7,
                label=f"{name}, tie lines")
    ax.set_xlabel("1000/T (1/K)")
    ax.set_ylabel("ln k")
    ax.legend(loc="upper left", frameon=False)
    core.save(fig, "fig4_vant_hoff")
    plt.close(fig)


if __name__ == "__main__":
    main()
