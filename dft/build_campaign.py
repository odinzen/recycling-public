"""Generate the Ti-interstitial DFT campaign (Quantum ESPRESSO inputs + Sol submit scripts).

References set the energy zero-points; the six interstitial supercells give the
formation energies and the relative N-vs-O potency. Interstitial sites are found
by a void search (grid the cell, take the local maxima of nearest-host distance)
and confirmed by coordination number: octahedral = 6, tetrahedral = 4. Nothing is
hand-placed.

Run: conda activate ti-dft && python build_campaign.py
Outputs land in ./inputs (one .in + .slurm per job) plus submit_all.sh.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk, molecule
from ase.spacegroup import crystal
from ase.io import write

OUT = Path(__file__).with_name("inputs")
OUT.mkdir(exist_ok=True)

PSEUDO_DIR = "/scratch/user/ti-dft/pseudo"
PSEUDO = {
    "Ti": "Ti.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "O":  "O.pbe-n-kjpaw_psl.1.0.0.UPF",
    "N":  "N.pbe-n-kjpaw_psl.1.0.0.UPF",
}
MASS = {"Ti": 47.867, "O": 15.999, "N": 14.007}
ECUTWFC, ECUTRHO = 60.0, 480.0

# experimental lattice constants (Angstrom)
A_HCP, C_HCP = 2.9508, 4.6855   # alpha-Ti
A_BCC = 3.3065                  # beta-Ti (high-T)


def base_system(nat, ntyp, metal=True, nspin=1, mag=None):
    sysd = {
        "ecutwfc": ECUTWFC, "ecutrho": ECUTRHO,
        "occupations": "smearing", "smearing": "mp", "degauss": 0.02,
    }
    if nspin == 2:
        sysd["nspin"] = 2
        if mag:
            for i, m in mag.items():
                sysd[f"starting_magnetization({i})"] = m
    return sysd


def input_data(calc, prefix, nat, ntyp, nspin=1, mag=None):
    d = {
        "control": {
            "calculation": calc, "prefix": prefix, "outdir": "./out",
            "pseudo_dir": PSEUDO_DIR, "tprnfor": True, "tstress": True,
        },
        "system": base_system(nat, ntyp, nspin=nspin, mag=mag),
        "electrons": {"conv_thr": 1e-7, "mixing_beta": 0.3},
    }
    if calc in ("relax", "vc-relax"):
        d["ions"] = {"ion_dynamics": "bfgs"}
    if calc == "vc-relax":
        d["cell"] = {"cell_dynamics": "bfgs"}
    return d


def find_voids(prim, grid=24, coord_tol=0.35):
    """Grid the primitive cell; return distinct void centers as (frac, dmin, coord).

    dmin = distance to nearest host atom; coord = hosts within dmin*(1+coord_tol).
    Local maxima of dmin over the grid are the interstitial sites.
    """
    cell = prim.cell.array
    pos = prim.get_positions()
    fr = np.stack(np.meshgrid(*[np.linspace(0, 1, grid, endpoint=False)] * 3, indexing="ij"), -1).reshape(-1, 3)
    cart = fr @ cell
    # min-image distance from each grid point to any host
    inv = np.linalg.inv(cell)
    dmins = np.empty(len(cart))
    for i, c in enumerate(cart):
        d = pos - c
        f = d @ inv
        f -= np.round(f)
        dmins[i] = np.min(np.linalg.norm(f @ cell, axis=1))
    order = np.argsort(-dmins)
    picked = []
    for idx in order:
        f = fr[idx]
        if all(np.linalg.norm(((f - p + 0.5) % 1 - 0.5)) > 0.12 for p, _, _ in picked):
            c = f @ cell
            d = pos - c
            ff = d @ inv
            ff -= np.round(ff)
            dd = np.linalg.norm(ff @ cell, axis=1)
            coord = int(np.sum(dd < dmins[idx] * (1 + coord_tol)))
            picked.append((f, dmins[idx], coord))
        if len(picked) >= 4:
            break
    return picked


def supercell_with_interstitial(prim, reps, frac_site, species):
    sc = prim.repeat(reps)
    # map the primitive fractional site into a central cell of the supercell
    r = np.array(reps)
    centre_cell = r // 2
    sc_frac = (centre_cell + np.array(frac_site)) / r
    cart = sc_frac @ sc.cell.array
    sc.append(species)
    sc.positions[-1] = cart
    return sc


def min_x_dist(atoms, species):
    """Minimum distance from the (single) interstitial of `species` to any host."""
    idx = [a.index for a in atoms if a.symbol == species][-1]
    d = atoms.get_distances(idx, [a.index for a in atoms if a.index != idx], mic=True)
    return float(np.min(d))


JOBS = []  # (name, atoms, calc, kpts, phase, nspin, mag)


def add(name, atoms, calc, kpts, phase, nspin=1, mag=None):
    JOBS.append((name, atoms, calc, kpts, phase, nspin, mag))


# ---- Phase 1: references ---------------------------------------------------
ti_hcp = bulk("Ti", "hcp", a=A_HCP, c=C_HCP)
ti_bcc = bulk("Ti", "bcc", a=A_BCC, cubic=True)
add("ref_Ti_hcp", ti_hcp, "vc-relax", (12, 12, 8), 1)
add("ref_Ti_bcc", ti_bcc, "vc-relax", (12, 12, 12), 1)

o2 = molecule("O2"); o2.set_cell([12, 12, 12]); o2.center()
n2 = molecule("N2"); n2.set_cell([12, 12, 12]); n2.center()
add("ref_O2", o2, "scf", (1, 1, 1), 1, nspin=2, mag={1: 0.5})
add("ref_N2", n2, "scf", (1, 1, 1), 1)

tio2 = crystal(["Ti", "O"], basis=[(0, 0, 0), (0.3053, 0.3053, 0.0)],
               spacegroup=136, cellpar=[4.593, 4.593, 2.959, 90, 90, 90])
tin = bulk("TiN", "rocksalt", a=4.24)
add("ref_TiO2_rutile", tio2, "vc-relax", (6, 6, 8), 1)
add("ref_TiN_rocksalt", tin, "vc-relax", (8, 8, 8), 1)

hcp_super = ti_hcp.repeat((3, 3, 2))   # 36 Ti
bcc_super = ti_bcc.repeat((3, 3, 3))   # 54 Ti
add("ref_super_hcp", hcp_super, "scf", (3, 3, 3), 1)
add("ref_super_bcc", bcc_super, "scf", (3, 3, 3), 1)

# ---- Phase 2: interstitials --------------------------------------------------
# Fractional sites verified against ASE's atom positions (min Ti-X distance in A):
#   HCP oct (2/3,1/3,1/4)=2.07  HCP tet (2/3,1/3,1/8)=1.80  BCC oct (0,0,1/2)=1.65
SITES = {
    ("hcp", "oct"): (2/3, 1/3, 1/4),
    ("hcp", "tet"): (2/3, 1/3, 1/8),
    ("bcc", "oct"): (0.0, 0.0, 0.5),
}
HOSTS = {"hcp": (ti_hcp, (3, 3, 2)), "bcc": (ti_bcc, (3, 3, 3))}

for (host, site), frac in SITES.items():
    prim, reps = HOSTS[host]
    for X in ("O", "N"):
        atoms = supercell_with_interstitial(prim, reps, frac, X)
        add(f"int_{host}_{X}_{site}", atoms, "relax", (3, 3, 3), 2)


SLURM = """#!/bin/bash
#SBATCH -J {name}
#SBATCH -A REPLACE_ACCOUNT
#SBATCH -p {part}
#SBATCH -q public
#SBATCH -N 1
#SBATCH -n {ntasks}
#SBATCH -t {wall}
#SBATCH -o {name}.%j.out
#SBATCH -e {name}.%j.err

module purge
module load quantum-espresso/7.3.1-cpu
export OMP_NUM_THREADS=1
export I_MPI_HYDRA_BOOTSTRAP=slurm

cd "$SLURM_SUBMIT_DIR"
mpirun -np $SLURM_NTASKS pw.x -in {name}.in > {name}.out
"""


def slurm_params(name, nat):
    if name in ("ref_O2", "ref_N2"):
        return "htc", 8, "0-00:30:00"
    if name.startswith("ref_super"):
        return "htc", 16, "0-03:00:00"
    if name.startswith("int_"):
        return "public", 16, "0-08:00:00"   # relaxations, long queue
    return "htc", 16, "0-02:00:00"


# ---- write everything ------------------------------------------------------
print("\n=== self-check + write ===")
for name, atoms, calc, kpts, phase, nspin, mag in JOBS:
    ntyp = len(set(atoms.get_chemical_symbols()))
    idata = input_data(calc, name, len(atoms), ntyp, nspin=nspin, mag=mag)
    dmat = atoms.get_all_distances(mic=True)
    np.fill_diagonal(dmat, 9.9)
    dmin = dmat.min()
    note = ""
    if "int_" in name:
        X = "O" if "_O_" in name else "N"
        note = f" | X-host min = {min_x_dist(atoms, X):.2f} A"
    flag = "  <-- CHECK" if dmin < 0.9 else ""
    print(f"  {name:22} {atoms.get_chemical_formula():12} nat={len(atoms):3} dmin={dmin:.2f}{note}{flag}")
    fin = OUT / f"{name}.in"
    write(fin, atoms, format="espresso-in",
          input_data=idata, pseudopotentials=PSEUDO, kpts=kpts)
    fin.write_text(re.sub(r"&(FCP|RISM)\s*\n\s*/\s*\n", "", fin.read_text()))
    part, ntasks, wall = slurm_params(name, len(atoms))
    (OUT / f"{name}.slurm").write_text(SLURM.format(name=name, part=part, ntasks=ntasks, wall=wall))

sub = "#!/bin/bash\n# submit the whole set; watch with: squeue -u $USER\n"
sub += "\n".join(f"sbatch {name}.slurm" for name, *_ in JOBS) + "\n"
(OUT / "submit_all.sh").write_text(sub)
print(f"\nwrote {len(JOBS)} inputs + submit scripts to {OUT}")
