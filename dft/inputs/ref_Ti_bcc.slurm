#!/bin/bash
#SBATCH -J ref_Ti_bcc
#SBATCH -A grp_klilova
#SBATCH -p htc
#SBATCH -q public
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 0-02:00:00
#SBATCH -o ref_Ti_bcc.%j.out
#SBATCH -e ref_Ti_bcc.%j.err

module purge
module load quantum-espresso/7.3.1-cpu
export OMP_NUM_THREADS=1
export I_MPI_HYDRA_BOOTSTRAP=slurm

cd "$SLURM_SUBMIT_DIR"
mpirun -np $SLURM_NTASKS pw.x -in ref_Ti_bcc.in > ref_Ti_bcc.out
