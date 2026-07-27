#!/bin/bash
#SBATCH -J ref_N2
#SBATCH -A grp_klilova
#SBATCH -p htc
#SBATCH -q public
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 0-00:30:00
#SBATCH -o ref_N2.%j.out
#SBATCH -e ref_N2.%j.err

module purge
module load quantum-espresso/7.3.1-cpu
export OMP_NUM_THREADS=1
export I_MPI_HYDRA_BOOTSTRAP=slurm

cd "$SLURM_SUBMIT_DIR"
mpirun -np $SLURM_NTASKS pw.x -in ref_N2.in > ref_N2.out
