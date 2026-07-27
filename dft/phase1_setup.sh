#!/bin/bash
# Phase-1 DFT reference jobs for the recycled-Ti study.
# Paste this whole block into your Sol terminal (it cd's itself).
mkdir -p /scratch/mebusta3/ti-dft/inputs
cd /scratch/mebusta3/ti-dft/inputs

cat > ref_N2.in <<'QEEOF'
&CONTROL
   calculation      = 'scf'
   tstress          = .true.
   tprnfor          = .true.
   outdir           = './out'
   prefix           = 'ref_N2'
   pseudo_dir       = '/scratch/mebusta3/ti-dft/pseudo'
/
&SYSTEM
   ecutwfc          = 60.0
   ecutrho          = 480.0
   occupations      = 'smearing'
   degauss          = 0.02
   smearing         = 'mp'
   ntyp             = 1
   nat              = 2
   ibrav            = 0
/
&ELECTRONS
   conv_thr         = 1e-07
   mixing_beta      = 0.3
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
N 14.007 N.pbe-n-kjpaw_psl.1.0.0.UPF

K_POINTS automatic
1 1 1  0 0 0

CELL_PARAMETERS angstrom
12.00000000000000 0.00000000000000 0.00000000000000
0.00000000000000 12.00000000000000 0.00000000000000
0.00000000000000 0.00000000000000 12.00000000000000

ATOMIC_POSITIONS angstrom
N 6.0000000000 6.0000000000 6.5649900000  
N 6.0000000000 6.0000000000 5.4350100000  
QEEOF

cat > ref_O2.in <<'QEEOF'
&CONTROL
   calculation      = 'scf'
   tstress          = .true.
   tprnfor          = .true.
   outdir           = './out'
   prefix           = 'ref_O2'
   pseudo_dir       = '/scratch/mebusta3/ti-dft/pseudo'
/
&SYSTEM
   starting_magnetization(1) = 1.0
   ecutwfc          = 60.0
   ecutrho          = 480.0
   occupations      = 'smearing'
   degauss          = 0.02
   smearing         = 'mp'
   nspin            = 2
   ntyp             = 1
   nat              = 2
   ibrav            = 0
/
&ELECTRONS
   conv_thr         = 1e-07
   mixing_beta      = 0.3
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
O  15.999 O.pbe-n-kjpaw_psl.1.0.0.UPF

K_POINTS automatic
1 1 1  0 0 0

CELL_PARAMETERS angstrom
12.00000000000000 0.00000000000000 0.00000000000000
0.00000000000000 12.00000000000000 0.00000000000000
0.00000000000000 0.00000000000000 12.00000000000000

ATOMIC_POSITIONS angstrom
O  6.0000000000 6.0000000000 6.6229780000  
O  6.0000000000 6.0000000000 5.3770220000  
QEEOF

cat > ref_super_bcc.in <<'QEEOF'
&CONTROL
   calculation      = 'scf'
   tstress          = .true.
   tprnfor          = .true.
   outdir           = './out'
   prefix           = 'ref_super_bcc'
   pseudo_dir       = '/scratch/mebusta3/ti-dft/pseudo'
/
&SYSTEM
   ecutwfc          = 60.0
   ecutrho          = 480.0
   occupations      = 'smearing'
   degauss          = 0.02
   smearing         = 'mp'
   ntyp             = 1
   nat              = 54
   ibrav            = 0
/
&ELECTRONS
   conv_thr         = 1e-07
   mixing_beta      = 0.3
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
Ti 47.867 Ti.pbe-spn-kjpaw_psl.1.0.0.UPF

K_POINTS automatic
3 3 3  0 0 0

CELL_PARAMETERS angstrom
9.91950000000000 0.00000000000000 0.00000000000000
0.00000000000000 9.91950000000000 0.00000000000000
0.00000000000000 0.00000000000000 9.91950000000000

ATOMIC_POSITIONS angstrom
Ti 0.0000000000 0.0000000000 0.0000000000  
Ti 1.6532500000 1.6532500000 1.6532500000  
Ti 0.0000000000 0.0000000000 3.3065000000  
Ti 1.6532500000 1.6532500000 4.9597500000  
Ti 0.0000000000 0.0000000000 6.6130000000  
Ti 1.6532500000 1.6532500000 8.2662500000  
Ti 0.0000000000 3.3065000000 0.0000000000  
Ti 1.6532500000 4.9597500000 1.6532500000  
Ti 0.0000000000 3.3065000000 3.3065000000  
Ti 1.6532500000 4.9597500000 4.9597500000  
Ti 0.0000000000 3.3065000000 6.6130000000  
Ti 1.6532500000 4.9597500000 8.2662500000  
Ti 0.0000000000 6.6130000000 0.0000000000  
Ti 1.6532500000 8.2662500000 1.6532500000  
Ti 0.0000000000 6.6130000000 3.3065000000  
Ti 1.6532500000 8.2662500000 4.9597500000  
Ti 0.0000000000 6.6130000000 6.6130000000  
Ti 1.6532500000 8.2662500000 8.2662500000  
Ti 3.3065000000 0.0000000000 0.0000000000  
Ti 4.9597500000 1.6532500000 1.6532500000  
Ti 3.3065000000 0.0000000000 3.3065000000  
Ti 4.9597500000 1.6532500000 4.9597500000  
Ti 3.3065000000 0.0000000000 6.6130000000  
Ti 4.9597500000 1.6532500000 8.2662500000  
Ti 3.3065000000 3.3065000000 0.0000000000  
Ti 4.9597500000 4.9597500000 1.6532500000  
Ti 3.3065000000 3.3065000000 3.3065000000  
Ti 4.9597500000 4.9597500000 4.9597500000  
Ti 3.3065000000 3.3065000000 6.6130000000  
Ti 4.9597500000 4.9597500000 8.2662500000  
Ti 3.3065000000 6.6130000000 0.0000000000  
Ti 4.9597500000 8.2662500000 1.6532500000  
Ti 3.3065000000 6.6130000000 3.3065000000  
Ti 4.9597500000 8.2662500000 4.9597500000  
Ti 3.3065000000 6.6130000000 6.6130000000  
Ti 4.9597500000 8.2662500000 8.2662500000  
Ti 6.6130000000 0.0000000000 0.0000000000  
Ti 8.2662500000 1.6532500000 1.6532500000  
Ti 6.6130000000 0.0000000000 3.3065000000  
Ti 8.2662500000 1.6532500000 4.9597500000  
Ti 6.6130000000 0.0000000000 6.6130000000  
Ti 8.2662500000 1.6532500000 8.2662500000  
Ti 6.6130000000 3.3065000000 0.0000000000  
Ti 8.2662500000 4.9597500000 1.6532500000  
Ti 6.6130000000 3.3065000000 3.3065000000  
Ti 8.2662500000 4.9597500000 4.9597500000  
Ti 6.6130000000 3.3065000000 6.6130000000  
Ti 8.2662500000 4.9597500000 8.2662500000  
Ti 6.6130000000 6.6130000000 0.0000000000  
Ti 8.2662500000 8.2662500000 1.6532500000  
Ti 6.6130000000 6.6130000000 3.3065000000  
Ti 8.2662500000 8.2662500000 4.9597500000  
Ti 6.6130000000 6.6130000000 6.6130000000  
Ti 8.2662500000 8.2662500000 8.2662500000  
QEEOF

cat > ref_super_hcp.in <<'QEEOF'
&CONTROL
   calculation      = 'scf'
   tstress          = .true.
   tprnfor          = .true.
   outdir           = './out'
   prefix           = 'ref_super_hcp'
   pseudo_dir       = '/scratch/mebusta3/ti-dft/pseudo'
/
&SYSTEM
   ecutwfc          = 60.0
   ecutrho          = 480.0
   occupations      = 'smearing'
   degauss          = 0.02
   smearing         = 'mp'
   ntyp             = 1
   nat              = 36
   ibrav            = 0
/
&ELECTRONS
   conv_thr         = 1e-07
   mixing_beta      = 0.3
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
Ti 47.867 Ti.pbe-spn-kjpaw_psl.1.0.0.UPF

K_POINTS automatic
3 3 3  0 0 0

CELL_PARAMETERS angstrom
8.85240000000000 0.00000000000000 0.00000000000000
-4.42620000000000 7.66640328446136 0.00000000000000
0.00000000000000 0.00000000000000 9.37100000000000

ATOMIC_POSITIONS angstrom
Ti 0.0000000000 0.0000000000 0.0000000000  
Ti 0.0000000000 1.7036451743 2.3427500000  
Ti 0.0000000000 0.0000000000 4.6855000000  
Ti 0.0000000000 1.7036451743 7.0282500000  
Ti -1.4754000000 2.5554677615 0.0000000000  
Ti -1.4754000000 4.2591129358 2.3427500000  
Ti -1.4754000000 2.5554677615 4.6855000000  
Ti -1.4754000000 4.2591129358 7.0282500000  
Ti -2.9508000000 5.1109355230 0.0000000000  
Ti -2.9508000000 6.8145806973 2.3427500000  
Ti -2.9508000000 5.1109355230 4.6855000000  
Ti -2.9508000000 6.8145806973 7.0282500000  
Ti 2.9508000000 0.0000000000 0.0000000000  
Ti 2.9508000000 1.7036451743 2.3427500000  
Ti 2.9508000000 0.0000000000 4.6855000000  
Ti 2.9508000000 1.7036451743 7.0282500000  
Ti 1.4754000000 2.5554677615 0.0000000000  
Ti 1.4754000000 4.2591129358 2.3427500000  
Ti 1.4754000000 2.5554677615 4.6855000000  
Ti 1.4754000000 4.2591129358 7.0282500000  
Ti 0.0000000000 5.1109355230 0.0000000000  
Ti 0.0000000000 6.8145806973 2.3427500000  
Ti 0.0000000000 5.1109355230 4.6855000000  
Ti 0.0000000000 6.8145806973 7.0282500000  
Ti 5.9016000000 0.0000000000 0.0000000000  
Ti 5.9016000000 1.7036451743 2.3427500000  
Ti 5.9016000000 0.0000000000 4.6855000000  
Ti 5.9016000000 1.7036451743 7.0282500000  
Ti 4.4262000000 2.5554677615 0.0000000000  
Ti 4.4262000000 4.2591129358 2.3427500000  
Ti 4.4262000000 2.5554677615 4.6855000000  
Ti 4.4262000000 4.2591129358 7.0282500000  
Ti 2.9508000000 5.1109355230 0.0000000000  
Ti 2.9508000000 6.8145806973 2.3427500000  
Ti 2.9508000000 5.1109355230 4.6855000000  
Ti 2.9508000000 6.8145806973 7.0282500000  
QEEOF

cat > ref_Ti_bcc.in <<'QEEOF'
&CONTROL
   calculation      = 'vc-relax'
   tstress          = .true.
   tprnfor          = .true.
   outdir           = './out'
   prefix           = 'ref_Ti_bcc'
   pseudo_dir       = '/scratch/mebusta3/ti-dft/pseudo'
/
&SYSTEM
   ecutwfc          = 60.0
   ecutrho          = 480.0
   occupations      = 'smearing'
   degauss          = 0.02
   smearing         = 'mp'
   ntyp             = 1
   nat              = 2
   ibrav            = 0
/
&ELECTRONS
   conv_thr         = 1e-07
   mixing_beta      = 0.3
/
&IONS
   ion_dynamics     = 'bfgs'
/
&CELL
   cell_dynamics    = 'bfgs'
/
ATOMIC_SPECIES
Ti 47.867 Ti.pbe-spn-kjpaw_psl.1.0.0.UPF

K_POINTS automatic
12 12 12  0 0 0

CELL_PARAMETERS angstrom
3.30650000000000 0.00000000000000 0.00000000000000
0.00000000000000 3.30650000000000 0.00000000000000
0.00000000000000 0.00000000000000 3.30650000000000

ATOMIC_POSITIONS angstrom
Ti 0.0000000000 0.0000000000 0.0000000000  
Ti 1.6532500000 1.6532500000 1.6532500000  
QEEOF

cat > ref_Ti_hcp.in <<'QEEOF'
&CONTROL
   calculation      = 'vc-relax'
   tstress          = .true.
   tprnfor          = .true.
   outdir           = './out'
   prefix           = 'ref_Ti_hcp'
   pseudo_dir       = '/scratch/mebusta3/ti-dft/pseudo'
/
&SYSTEM
   ecutwfc          = 60.0
   ecutrho          = 480.0
   occupations      = 'smearing'
   degauss          = 0.02
   smearing         = 'mp'
   ntyp             = 1
   nat              = 2
   ibrav            = 0
/
&ELECTRONS
   conv_thr         = 1e-07
   mixing_beta      = 0.3
/
&IONS
   ion_dynamics     = 'bfgs'
/
&CELL
   cell_dynamics    = 'bfgs'
/
ATOMIC_SPECIES
Ti 47.867 Ti.pbe-spn-kjpaw_psl.1.0.0.UPF

K_POINTS automatic
12 12 8  0 0 0

CELL_PARAMETERS angstrom
2.95080000000000 0.00000000000000 0.00000000000000
-1.47540000000000 2.55546776148712 0.00000000000000
0.00000000000000 0.00000000000000 4.68550000000000

ATOMIC_POSITIONS angstrom
Ti 0.0000000000 0.0000000000 0.0000000000  
Ti 0.0000000000 1.7036451743 2.3427500000  
QEEOF

cat > ref_TiN_rocksalt.in <<'QEEOF'
&CONTROL
   calculation      = 'vc-relax'
   tstress          = .true.
   tprnfor          = .true.
   outdir           = './out'
   prefix           = 'ref_TiN_rocksalt'
   pseudo_dir       = '/scratch/mebusta3/ti-dft/pseudo'
/
&SYSTEM
   ecutwfc          = 60.0
   ecutrho          = 480.0
   occupations      = 'smearing'
   degauss          = 0.02
   smearing         = 'mp'
   ntyp             = 2
   nat              = 2
   ibrav            = 0
/
&ELECTRONS
   conv_thr         = 1e-07
   mixing_beta      = 0.3
/
&IONS
   ion_dynamics     = 'bfgs'
/
&CELL
   cell_dynamics    = 'bfgs'
/
ATOMIC_SPECIES
Ti 47.867 Ti.pbe-spn-kjpaw_psl.1.0.0.UPF
N 14.007 N.pbe-n-kjpaw_psl.1.0.0.UPF

K_POINTS automatic
8 8 8  0 0 0

CELL_PARAMETERS angstrom
0.00000000000000 2.12000000000000 2.12000000000000
2.12000000000000 0.00000000000000 2.12000000000000
2.12000000000000 2.12000000000000 0.00000000000000

ATOMIC_POSITIONS angstrom
Ti 0.0000000000 0.0000000000 0.0000000000  
N 2.1200000000 0.0000000000 0.0000000000  
QEEOF

cat > ref_TiO2_rutile.in <<'QEEOF'
&CONTROL
   calculation      = 'vc-relax'
   tstress          = .true.
   tprnfor          = .true.
   outdir           = './out'
   prefix           = 'ref_TiO2_rutile'
   pseudo_dir       = '/scratch/mebusta3/ti-dft/pseudo'
/
&SYSTEM
   ecutwfc          = 60.0
   ecutrho          = 480.0
   occupations      = 'smearing'
   degauss          = 0.02
   smearing         = 'mp'
   ntyp             = 2
   nat              = 6
   ibrav            = 0
/
&ELECTRONS
   conv_thr         = 1e-07
   mixing_beta      = 0.3
/
&IONS
   ion_dynamics     = 'bfgs'
/
&CELL
   cell_dynamics    = 'bfgs'
/
ATOMIC_SPECIES
Ti 47.867 Ti.pbe-spn-kjpaw_psl.1.0.0.UPF
O 15.999 O.pbe-n-kjpaw_psl.1.0.0.UPF

K_POINTS automatic
6 6 8  0 0 0

CELL_PARAMETERS angstrom
4.59300000000000 0.00000000000000 0.00000000000000
0.00000000000000 4.59300000000000 0.00000000000000
0.00000000000000 0.00000000000000 2.95900000000000

ATOMIC_POSITIONS angstrom
Ti 0.0000000000 0.0000000000 0.0000000000  
Ti 2.2965000000 2.2965000000 1.4795000000  
O 1.4022429000 1.4022429000 0.0000000000  
O 3.1907571000 3.1907571000 0.0000000000  
O 0.8942571000 3.6987429000 1.4795000000  
O 3.6987429000 0.8942571000 1.4795000000  
QEEOF

for n in ref_Ti_hcp ref_Ti_bcc ref_O2 ref_N2 ref_TiO2_rutile ref_TiN_rocksalt ref_super_hcp ref_super_bcc; do
  nt=16; wt="0-02:00:00"
  case $n in ref_O2|ref_N2) nt=8; wt="0-00:30:00";; ref_super_*) wt="0-03:00:00";; esac
  cat > $n.slurm <<SLEOF
#!/bin/bash
#SBATCH -J $n
#SBATCH -A grp_klilova
#SBATCH -p htc
#SBATCH -q public
#SBATCH -N 1
#SBATCH -n $nt
#SBATCH -t $wt
#SBATCH -o $n.%j.out
#SBATCH -e $n.%j.err

module purge
module load quantum-espresso/7.3.1-cpu
export OMP_NUM_THREADS=1
export I_MPI_HYDRA_BOOTSTRAP=slurm
cd "\$SLURM_SUBMIT_DIR"
mpirun -np \$SLURM_NTASKS pw.x -in $n.in > $n.out
SLEOF
done
echo
echo "Ready: $(ls *.in 2>/dev/null | wc -l) inputs, $(ls *.slurm 2>/dev/null | wc -l) submit scripts in $(pwd)"
