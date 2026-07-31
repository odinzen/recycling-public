#!/bin/bash
#SBATCH -J int_hcp_O_oct
#SBATCH -A REPLACE_ACCOUNT
#SBATCH -p public
#SBATCH -q public
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 0-08:00:00
#SBATCH -o int_hcp_O_oct.%j.out
#SBATCH -e int_hcp_O_oct.%j.err

module purge
module load quantum-espresso/7.3.1-cpu
export OMP_NUM_THREADS=1
export I_MPI_HYDRA_BOOTSTRAP=slurm

cd "$SLURM_SUBMIT_DIR"
mpirun -np $SLURM_NTASKS pw.x -in int_hcp_O_oct.in > int_hcp_O_oct.out
