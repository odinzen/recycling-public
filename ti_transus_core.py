"""Calculations shared by the Ti-O and Ti-N scripts: transus, tie lines, and the dilute alpha-beta transfer terms.

The transus is found by bisection on temperature (to 0.01 K) rather than read off a temperature step, since a 5 K step
had put pure Ti at 1152.5 K and quantized every slope.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from pycalphad import Database, calculate, equilibrium, variables as v

DATA = Path(__file__).resolve().parent / "data"
TDB = {"O": DATA / "ti_o_lindwall_metallic_v1.tdb", "N": DATA / "ti_n_ohtani_metallic_v2.tdb"}
MASS = {"O": 15.999, "N": 14.007, "TI": 47.867}
R = 8.314462618
PHASES = ["HCP_A3", "BCC_A2"]
SITES = {"HCP_A3": 0.5, "BCC_A2": 3.0}  # interstitial sites per Ti in each reduced database


def load(el: str) -> Database:
    return Database(str(TDB[el]))


def wt_pct(x: float, el: str) -> float:
    return 100.0 * x * MASS[el] / (x * MASS[el] + (1 - x) * MASS["TI"])


def _equilibrium(db, el, x, T):
    if x < 1e-12:
        return equilibrium(db, ["TI", "VA"], PHASES, {v.T: T, v.P: 101325, v.N: 1})
    return equilibrium(db, ["TI", el, "VA"], PHASES, {v.X(el): x, v.T: T, v.P: 101325, v.N: 1})


def hcp_fraction(db, el, x, T) -> float:
    eq = _equilibrium(db, el, x, T)
    ph, n = eq.Phase.values.squeeze(), eq.NP.values.squeeze()
    return float(sum(nn for p, nn in zip(ph, n) if p == "HCP_A3" and nn == nn))


def transus(db, el, x, lo=1100.0, hi=1450.0, tol=0.01) -> float:
    """Beta transus (K) at mole fraction x: the temperature at which the last alpha disappears on heating."""
    if hcp_fraction(db, el, x, lo) < 1e-6 or hcp_fraction(db, el, x, hi) > 1e-6:
        raise ValueError(f"transus for {el} x={x} is outside {lo} to {hi} K")
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if hcp_fraction(db, el, x, mid) > 1e-6 else (lo, mid)
    return 0.5 * (lo + hi)


def tie_line(db, el, T, x_max=0.05, n=200):
    """Mole fractions of el in the coexisting alpha and beta at T (K), or None outside the two-phase field."""
    xs = np.linspace(2e-4, x_max, n)
    eq = equilibrium(db, ["TI", el, "VA"], PHASES, {v.X(el): xs, v.T: T, v.P: 101325, v.N: 1})
    ph = eq.Phase.values.squeeze()
    X = eq.X.sel(component=el).values.squeeze()
    NP = eq.NP.values.squeeze()
    for i in range(n):
        present = {p: (nn, xx) for p, nn, xx in zip(ph[i], NP[i], X[i]) if p in PHASES}
        if len(present) == 2 and min(present["HCP_A3"][0], present["BCC_A2"][0]) > 0.05:
            return float(present["HCP_A3"][1]), float(present["BCC_A2"][1])
    return None


def _formula_gibbs(db, el, phase, y, T):
    """Gibbs energy per formula unit Ti(el_y Va_1-y)_a, J."""
    r = calculate(db, ["TI", el, "VA"], phase, T=T, P=101325, N=1, points=np.array([[1.0, y, 1 - y]]), output="GM")
    return float(r.GM.values.squeeze()) * (1 + SITES[phase] * y)


def partial_gibbs_dilute(db, el, phase, T) -> float:
    """Dilute-limit partial Gibbs energy of el in a phase, less its configurational term, J/mol el.

    For (Ti)1(X,Va)a, mu_X = [G(Ti:X) - G(Ti:Va)]/a + L0/a + RT ln(y) as y -> 0; this returns the first two terms.
    L0 is read back from the database by evaluating the phase at y = 1/2.
    """
    a = SITES[phase]
    gx, gv, gh = (_formula_gibbs(db, el, phase, y, T) for y in (1.0, 0.0, 0.5))
    L0 = 4 * (gh - 0.5 * (gx + gv) - a * R * T * np.log(0.5))
    return (gx - gv) / a + L0 / a


def transfer_gibbs(db, el, T) -> float:
    """Dilute alpha-minus-beta partial Gibbs energy of el, configurational term excluded, J/mol el."""
    return partial_gibbs_dilute(db, el, "HCP_A3", T) - partial_gibbs_dilute(db, el, "BCC_A2", T)


def dilute_partition(db, el, T) -> float:
    """Dilute partition coefficient k = (el/Ti in alpha)/(el/Ti in beta)."""
    return float(np.exp(np.log(SITES["HCP_A3"] / SITES["BCC_A2"]) - transfer_gibbs(db, el, T) / (R * T)))


def transfer_split(db, el, T=1180.0, dT=10.0) -> dict:
    """Enthalpy and entropy of moving one mole of el from beta to alpha at infinite dilution, at T.

    The configurational part, R ln(a_alpha / a_beta), is the same for every interstitial in these two databases.
    """
    g_lo, g_hi = transfer_gibbs(db, el, T - dT), transfer_gibbs(db, el, T + dT)
    s_nc = -(g_hi - g_lo) / (2 * dT)
    g = transfer_gibbs(db, el, T)
    s_conf = R * np.log(SITES["HCP_A3"] / SITES["BCC_A2"])
    return {"T": T, "dH": g + T * s_nc, "dS_nonconfig": s_nc, "dS_config": s_conf, "dS_total": s_nc + s_conf,
            "k_dilute": dilute_partition(db, el, T)}


def ti_transformation_enthalpy(db, T=1155.0) -> float:
    """H(bcc) - H(hcp) of pure Ti, J/mol, from the database's own lattice stabilities."""
    def h(phase):
        r = calculate(db, ["TI", "VA"], phase, T=T, P=101325, N=1, points=np.array([[1.0, 1.0]]), output="HM")
        return float(r.HM.values.squeeze())
    return h("BCC_A2") - h("HCP_A3")


def figure_axes(figsize=(5.5, 4.2), ncols=1):
    """Greyscale house style: black frame, no grid, fonts readable at single-column width."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.linewidth": 1.0,
                         "legend.fontsize": 10, "xtick.direction": "in", "ytick.direction": "in"})
    fig, axes = plt.subplots(1, ncols, figsize=figsize)
    for ax in np.atleast_1d(axes):
        ax.set_facecolor("white")
        for s in ax.spines.values():
            s.set_color("black")
    return plt, fig, axes


def save(fig, name):
    out = Path(__file__).resolve().parent / "figures" / name
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    return out
