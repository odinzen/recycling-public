#!/usr/bin/env python
"""Self-contained Ti-interstitial DFT input generator (pure stdlib, no ASE/numpy).

Paste-and-run on Sol: writes all 14 Quantum ESPRESSO inputs + SLURM scripts into
./inputs. Structures and interstitial sites match the locally verified set
(HCP oct/tet, BCC oct; Ti-X starts 2.07/1.80/1.65 A).
"""
import math, os

A_HCP, C_HCP = 2.9508, 4.6855
A_BCC = 3.3065
PDIR = "/scratch/mebusta3/ti-dft/pseudo"
PS = {"Ti": "Ti.pbe-spn-kjpaw_psl.1.0.0.UPF",
      "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
      "N": "N.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS = {"Ti": 47.867, "O": 15.999, "N": 14.007}
OUT = "inputs"
os.makedirs(OUT, exist_ok=True)


def fc(cell, f):  # fractional -> cartesian
    return [sum(f[k] * cell[k][i] for k in range(3)) for i in range(3)]


def hexcell(a, c):
    return [[a, 0, 0], [-a / 2, a * math.sqrt(3) / 2, 0], [0, 0, c]]


def cub(a):
    return [[a, 0, 0], [0, a, 0], [0, 0, a]]


def supercell(cell, basis, reps):
    nx, ny, nz = reps
    sc = [[cell[i][j] * (reps[i]) for j in range(3)] for i in range(3)]
    atoms = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for sym, f in basis:
                    ff = [(f[0] + i) / nx, (f[1] + j) / ny, (f[2] + k) / nz]
                    atoms.append((sym, fc(sc, ff)))
    return sc, atoms


def add_interstitial(sc, atoms, reps, frac, X):
    ctr = [r // 2 for r in reps]
    ff = [(ctr[m] + frac[m]) / reps[m] for m in range(3)]
    return atoms + [(X, fc(sc, ff))]


def mindist(cell, atoms):
    # min pairwise distance with min-image (orthorhombic-safe via 27 images)
    import itertools
    best = 9.9
    imgs = list(itertools.product((-1, 0, 1), repeat=3))
    for a in range(len(atoms)):
        for b in range(a + 1, len(atoms)):
            for im in imgs:
                sh = [sum(im[k] * cell[k][i] for k in range(3)) for i in range(3)]
                d = sum((atoms[a][1][i]-atoms[b][1][i]-sh[i])**2 for i in range(3))**0.5
                best = min(best, d)
    return best


def write_qe(name, cell, atoms, kpts, calc="scf", nspin=1, mag=None):
    syms = []
    for s, _ in atoms:
        if s not in syms:
            syms.append(s)
    L = ["&CONTROL",
         "   calculation = '%s'" % calc,
         "   tstress = .true.", "   tprnfor = .true.",
         "   outdir = './out'", "   prefix = '%s'" % name,
         "   pseudo_dir = '%s'" % PDIR, "/", "&SYSTEM",
         "   ibrav = 0", "   nat = %d" % len(atoms), "   ntyp = %d" % len(syms),
         "   ecutwfc = 60.0", "   ecutrho = 480.0",
         "   occupations = 'smearing'", "   smearing = 'mp'", "   degauss = 0.02"]
    if nspin == 2:
        L.append("   nspin = 2")
        L.append("   starting_magnetization(1) = %.2f" % (mag or 0.5))
    L += ["/", "&ELECTRONS", "   conv_thr = 1d-7", "   mixing_beta = 0.3", "/"]
    if calc in ("relax", "vc-relax"):
        L += ["&IONS", "   ion_dynamics = 'bfgs'", "/"]
    if calc == "vc-relax":
        L += ["&CELL", "   cell_dynamics = 'bfgs'", "/"]
    L.append("ATOMIC_SPECIES")
    for s in syms:
        L.append("   %-3s %8.3f  %s" % (s, MASS[s], PS[s]))
    L.append("CELL_PARAMETERS angstrom")
    for v in cell:
        L.append("   %18.12f %18.12f %18.12f" % (v[0], v[1], v[2]))
    L.append("ATOMIC_POSITIONS angstrom")
    for s, p in atoms:
        L.append("   %-3s %16.10f %16.10f %16.10f" % (s, p[0], p[1], p[2]))
    L.append("K_POINTS automatic")
    L.append("   %d %d %d 0 0 0" % kpts)
    open(os.path.join(OUT, name + ".in"), "w").write("\n".join(L) + "\n")


def write_slurm(name):
    part, nt, wall = ("public", 16, "0-08:00:00") if name.startswith("int_") else \
        ("htc", 8, "0-00:30:00") if name in ("ref_O2", "ref_N2") else \
        ("htc", 16, "0-03:00:00") if name.startswith("ref_super") else ("htc", 16, "0-02:00:00")
    s = ("#!/bin/bash\n#SBATCH -J %s\n#SBATCH -A grp_klilova\n#SBATCH -p %s\n"
         "#SBATCH -q public\n#SBATCH -N 1\n#SBATCH -n %d\n#SBATCH -t %s\n"
         "#SBATCH -o %s.%%j.out\n#SBATCH -e %s.%%j.err\n\n"
         "module purge\nmodule load quantum-espresso/7.3.1-cpu\n"
         "export OMP_NUM_THREADS=1\nexport I_MPI_HYDRA_BOOTSTRAP=slurm\n\n"
         'cd "$SLURM_SUBMIT_DIR"\nmpirun -np $SLURM_NTASKS pw.x -in %s.in > %s.out\n'
         % (name, part, nt, wall, name, name, name, name))
    open(os.path.join(OUT, name + ".slurm"), "w").write(s)


# ---- structures ----
HCP = (hexcell(A_HCP, C_HCP), [("Ti", (0, 0, 0)), ("Ti", (1 / 3, 2 / 3, 1 / 2))])
BCC = (cub(A_BCC), [("Ti", (0, 0, 0)), ("Ti", (0.5, 0.5, 0.5))])
JOBS = []

# references
JOBS.append(("ref_Ti_hcp", HCP[0], [(s, fc(HCP[0], f)) for s, f in HCP[1]], (12, 12, 8), "vc-relax", 1, None))
JOBS.append(("ref_Ti_bcc", BCC[0], [(s, fc(BCC[0], f)) for s, f in BCC[1]], (12, 12, 12), "vc-relax", 1, None))
box = cub(12.0)
JOBS.append(("ref_O2", box, [("O", [6, 6, 6 + 0.6115]), ("O", [6, 6, 6 - 0.6115])], (1, 1, 1), "scf", 2, 1.0))
JOBS.append(("ref_N2", box, [("N", [6, 6, 6 + 0.55]), ("N", [6, 6, 6 - 0.55])], (1, 1, 1), "scf", 1, None))
# rutile TiO2 (P4_2/mnm), u=0.3053
ar, cr, u = 4.593, 2.959, 0.3053
rut = cub(ar); rut[2][2] = cr
tio2 = [("Ti", (0, 0, 0)), ("Ti", (0.5, 0.5, 0.5)),
        ("O", (u, u, 0)), ("O", (-u, -u, 0)), ("O", (0.5 + u, 0.5 - u, 0.5)), ("O", (0.5 - u, 0.5 + u, 0.5))]
JOBS.append(("ref_TiO2_rutile", rut, [(s, fc(rut, f)) for s, f in tio2], (6, 6, 8), "vc-relax", 1, None))
# rocksalt TiN conventional 8-atom
an = 4.24; rs = cub(an)
tin = [("Ti", (0, 0, 0)), ("Ti", (0.5, 0.5, 0)), ("Ti", (0.5, 0, 0.5)), ("Ti", (0, 0.5, 0.5)),
       ("N", (0.5, 0, 0)), ("N", (0, 0.5, 0)), ("N", (0, 0, 0.5)), ("N", (0.5, 0.5, 0.5))]
JOBS.append(("ref_TiN_rocksalt", rs, [(s, fc(rs, f)) for s, f in tin], (6, 6, 6), "vc-relax", 1, None))
# pure supercells
sc_h, at_h = supercell(HCP[0], HCP[1], (3, 3, 2))
sc_b, at_b = supercell(BCC[0], BCC[1], (3, 3, 3))
JOBS.append(("ref_super_hcp", sc_h, at_h, (3, 3, 3), "scf", 1, None))
JOBS.append(("ref_super_bcc", sc_b, at_b, (3, 3, 3), "scf", 1, None))

# interstitials (verified fractional sites)
SITES = {("hcp", "oct"): (2 / 3, 1 / 3, 1 / 4), ("hcp", "tet"): (2 / 3, 1 / 3, 1 / 8),
         ("bcc", "oct"): (0.0, 0.0, 0.5)}
for (host, site), frac in SITES.items():
    (cell, basis), reps = (HCP, (3, 3, 2)) if host == "hcp" else (BCC, (3, 3, 3))
    sc, base = supercell(cell, basis, reps)
    for X in ("O", "N"):
        at = add_interstitial(sc, base, reps, frac, X)
        JOBS.append(("int_%s_%s_%s" % (host, X, site), sc, at, (3, 3, 3), "relax", 1, None))

print("name                    nat  dmin(A)")
for name, cell, atoms, kpts, calc, nspin, mag in JOBS:
    write_qe(name, cell, atoms, kpts, calc, nspin, mag)
    write_slurm(name)
    tag = "  X-host=%.2f" % mindist(cell, atoms) if name.startswith("int_") else ""
    print("  %-22s %3d  %.2f%s" % (name, len(atoms), mindist(cell, atoms), tag))
print("\nWrote %d inputs + SLURM scripts to ./%s" % (len(JOBS), OUT))
