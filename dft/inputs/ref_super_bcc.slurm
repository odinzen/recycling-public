#!/bin/bash
#SBATCH -J ref_super_bcc
#SBATCH -A grp_klilova
#SBATCH -p htc
#SBATCH -q public
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 0-03:00:00
#SBATCH -o ref_super_bcc.%j.out
#SBATCH -e ref_super_bcc.%j.err

module purge
module load quantum-espresso/7.3.1-cpu
export OMP_NUM_THREADS=1
export I_MPI_HYDRA_BOOTSTRAP=slurm

cd "$SLURM_SUBMIT_DIR"
mpirun -np $SLURM_NTASKS pw.x -in ref_super_bcc.in > ref_super_bcc.out
