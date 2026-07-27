#!/bin/bash
# submit the whole set; watch with: squeue -u $USER
sbatch ref_Ti_hcp.slurm
sbatch ref_Ti_bcc.slurm
sbatch ref_O2.slurm
sbatch ref_N2.slurm
sbatch ref_TiO2_rutile.slurm
sbatch ref_TiN_rocksalt.slurm
sbatch ref_super_hcp.slurm
sbatch ref_super_bcc.slurm
sbatch int_hcp_O_oct.slurm
sbatch int_hcp_N_oct.slurm
sbatch int_hcp_O_tet.slurm
sbatch int_hcp_N_tet.slurm
sbatch int_bcc_O_oct.slurm
sbatch int_bcc_N_oct.slurm
