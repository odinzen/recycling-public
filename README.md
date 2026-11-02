# Recycled aerospace titanium: thermodynamic and first-principles data

Reproducibility package (analysis code and density functional data) for a CALPHAD
and first-principles study of interstitial oxygen and nitrogen in titanium, motivated
by the qualification of recycled aerospace Ti. The beta transus of Ti-O and Ti-N is
calculated from reduced metallic CALPHAD databases. Nitrogen raises the transus about
1.6 times as much as oxygen per wt%, below the factor of two of the oxygen equivalent
O_eq = O + 2*N + 0.67*C, and a dilute-solution analysis of the same databases traces
the difference to the solution entropy of the two interstitials.

This repository holds the analysis code and the DFT data only. The accompanying
manuscript is in preparation and will be linked on publication.

## CALPHAD analysis (Python, pycalphad)
- `ti_transus_core.py`          shared transus (bisection), tie lines, dilute transfer terms
- `ti_validation.py`            pure-Ti and Ti-V transus (validation)
- `ti_oeq_coefficient_test.py`  oxygen and nitrogen transus, nitrogen:oxygen weighting
- `ti_partition.py`             alpha/beta tie lines and partition coefficients
- `ti_transfer_entropy.py`      enthalpy and entropy of moving O or N from beta to alpha
- `reproduce.sh`                runs the sequence above

Each script writes a `*_results.json` holding the numbers behind the paper's tables and figures.
`ti_interstitial_params.json` is the sourced parameter ledger; every value carries a citation.

## First-principles data (Quantum ESPRESSO)
`dft/` holds the input generators and the exact QE input and SLURM files for the
campaign: reference states (Ti hcp/bcc, O2, N2, TiO2, TiN, pure supercells) and
oxygen and nitrogen interstitials at octahedral and tetrahedral sites in HCP and
BCC titanium, plus the k-point convergence batch (`kconv_on_sol.py`). The
calculations were run on the Sol supercomputer at Arizona State University.

## Thermodynamic databases
The reduced Ti-O and Ti-N metallic databases are own-derived reconstructions of
published assessments. They are available from the corresponding author on
reasonable request and are not distributed here; the compute scripts expect them
under `data/`.

## Environment
    conda env create -f environment.yml
    conda activate ti_recycle_env
    bash reproduce.sh

## License
The code is under the MIT License (see `LICENSE`); the data and figures are under CC BY 4.0 (see `LICENSE-DATA`). Copyright 2026 Odinzen LLC. Cite the paper if you use either.
