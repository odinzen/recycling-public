#!/usr/bin/env python
"""k-point convergence batch for the alpha-vs-beta interstitial preference.

Reuses the RELAXED geometries from the finished int_*.out (single-point SCF, no
re-relaxation) plus freshly built pure supercells, at denser k-meshes. Writes to
./kconv. Run on Sol from /scratch/user/ti-dft/inputs (needs the *.out there).
"""
import math, os, re

A_HCP, C_HCP = 2.9508, 4.6855
A_BCC = 3.3065
PDIR = "/scratch/user/ti-dft/pseudo"
PS = {"Ti": "Ti.pbe-spn-kjpaw_psl.1.0.0.UPF",
      "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF", "N": "N.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS = {"Ti": 47.867, "O": 15.999, "N": 14.007}
OUT = "kconv"
os.makedirs(OUT, exist_ok=True)
KMESHES = [(4, 4, 4), (5, 5, 5), (6, 6, 6)]


def fc(cell, f):
    return [sum(f[k] * cell[k][i] for k in range(3)) for i in range(3)]


def hexcell(a, c):
    return [[a, 0, 0], [-a / 2, a * math.sqrt(3) / 2, 0], [0, 0, c]]


def cub(a):
    return [[a, 0, 0], [0, a, 0], [0, 0, a]]


def supercell(cell, basis, reps):
    nx, ny, nz = reps
    sc = [[cell[i][j] * reps[i] for j in range(3)] for i in range(3)]
    atoms = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for sym, f in basis:
                    atoms.append((sym, fc(sc, [(f[0] + i) / nx, (f[1] + j) / ny, (f[2] + k) / nz])))
    return sc, atoms


def parse_relaxed(outfile):
    """Return [(sym,[x,y,z]), ...] from the 'Begin final coordinates' block;
    fall back to the last ATOMIC_POSITIONS block. Angstrom."""
    txt = open(outfile).read()
    m = re.search(r"Begin final coordinates(.*?)End final coordinates", txt, re.S)
    blk = m.group(1) if m else txt
    aps = list(re.finditer(r"ATOMIC_POSITIONS[^\n]*\n(.*?)(?:\n\s*\n|\Z)", blk, re.S))
    if not aps:
        raise RuntimeError("no ATOMIC_POSITIONS found in %s" % outfile)
    atoms = []
    for line in aps[-1].group(1).splitlines():
        p = line.split()
        if len(p) >= 4 and p[0] in MASS:
            atoms.append((p[0], [float(p[1]), float(p[2]), float(p[3])]))
    return atoms


def write_scf(name, cell, atoms, kpts):
    syms = []
    for s, _ in atoms:
        if s not in syms:
            syms.append(s)
    nsp = 2 if "O" in syms and False else 1  # interstitial cells are non-magnetic
    L = ["&CONTROL", "   calculation = 'scf'", "   outdir = './out'",
         "   prefix = '%s'" % name, "   pseudo_dir = '%s'" % PDIR, "/",
         "&SYSTEM", "   ibrav = 0", "   nat = %d" % len(atoms), "   ntyp = %d" % len(syms),
         "   ecutwfc = 60.0", "   ecutrho = 480.0",
         "   occupations = 'smearing'", "   smearing = 'mp'", "   degauss = 0.02", "/",
         "&ELECTRONS", "   conv_thr = 1d-7", "   mixing_beta = 0.3", "/",
         "ATOMIC_SPECIES"]
    for s in syms:
        L.append("   %-3s %8.3f  %s" % (s, MASS[s], PS[s]))
    L.append("CELL_PARAMETERS angstrom")
    for v in cell:
        L.append("   %18.12f %18.12f %18.12f" % (v[0], v[1], v[2]))
    L.append("ATOMIC_POSITIONS angstrom")
    for s, p in atoms:
        L.append("   %-3s %16.10f %16.10f %16.10f" % (s, p[0], p[1], p[2]))
    L += ["K_POINTS automatic", "   %d %d %d 0 0 0" % kpts]
    open(os.path.join(OUT, name + ".in"), "w").write("\n".join(L) + "\n")


def write_slurm(name, nat):
    nt, wall = (32, "0-08:00:00") if nat > 40 else (16, "0-04:00:00")
    s = ("#!/bin/bash\n#SBATCH -J %s\n#SBATCH -A REPLACE_ACCOUNT\n#SBATCH -p public\n"
         "#SBATCH -q public\n#SBATCH -N 1\n#SBATCH -n %d\n#SBATCH -t %s\n"
         "#SBATCH -o %s.%%j.out\n#SBATCH -e %s.%%j.err\n\n"
         "module purge\nmodule load quantum-espresso/7.3.1-cpu\n"
         "export OMP_NUM_THREADS=1\nexport I_MPI_HYDRA_BOOTSTRAP=slurm\n\n"
         'cd "$SLURM_SUBMIT_DIR"\nmpirun -np $SLURM_NTASKS pw.x -in %s.in > %s.out\n'
         % (name, nt, wall, name, name, name, name))
    open(os.path.join(OUT, name + ".slurm"), "w").write(s)


# ---- structures ----
HCP = (hexcell(A_HCP, C_HCP), [("Ti", (0, 0, 0)), ("Ti", (1 / 3, 2 / 3, 1 / 2))])
BCC = (cub(A_BCC), [("Ti", (0, 0, 0)), ("Ti", (0.5, 0.5, 0.5))])
hcp_cell, hcp_atoms = supercell(HCP[0], HCP[1], (3, 3, 2))
bcc_cell, bcc_atoms = supercell(BCC[0], BCC[1], (3, 3, 3))

STRUCTS = {"hcp_pure": (hcp_cell, hcp_atoms), "bcc_pure": (bcc_cell, bcc_atoms)}
# relaxed interstitials from the finished runs (single-point SCF at their final geometry)
for tag, cell, src in [("hcp_O", hcp_cell, "int_hcp_O_oct.out"),
                       ("hcp_N", hcp_cell, "int_hcp_N_oct.out"),
                       ("bcc_O", bcc_cell, "int_bcc_O_oct.out"),
                       ("bcc_N", bcc_cell, "int_bcc_N_oct.out")]:
    STRUCTS[tag] = (cell, parse_relaxed(src))

print("structure               nat")
n = 0
for tag, (cell, atoms) in STRUCTS.items():
    print("  %-20s %3d" % (tag, len(atoms)))
    for kx, ky, kz in KMESHES:
        name = "%s_k%d" % (tag, kx)
        write_scf(name, cell, atoms, (kx, ky, kz))
        write_slurm(name, len(atoms))
        n += 1
print("\nWrote %d SCF inputs (+ slurm) to ./%s" % (n, OUT))
