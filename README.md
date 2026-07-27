# Recycled aerospace titanium: thermodynamic and first-principles data

Reproducibility package (analysis code and density functional data) for a CALPHAD
and first-principles study of interstitial oxygen and nitrogen in titanium, aimed
at qualifying fully recycled aerospace Ti (Ti-6Al-4V, Ti-900 class). The oxygen
equivalent O_eq = O + 2*N + 0.67*C is the controlling acceptance variable, and the
alpha/beta transus is computed from reduced metallic CALPHAD databases. The central
result is that nitrogen's alpha-stabilizing potency is about 1.5 times oxygen's,
below the factor of two the oxygen equivalent assigns.

This repository holds the analysis code and the DFT data only. The accompanying
manuscript is in preparation and will be linked on publication.

## CALPHAD analysis (Python, pycalphad)
- `ti_validation.py`            pure-Ti and Ti-V transus (validation)
- `ti_o_transus.py`             Ti-O transus shift with oxygen
- `ti_o_partition.py`           Ti-O alpha/beta partition and two-phase field
- `ti_oeq_coefficient_test.py`  nitrogen vs oxygen alpha-stabilizing potency
- `ti_transus_figures.py`       schematic figures
- `reproduce.sh`                runs the sequence above

Each script writes a `*_results.json` with a provenance header.
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
under `data/recycled_ti/`.

## Environment
    conda env create -f environment.yml
    conda activate ti_recycle_env
    bash reproduce.sh

## License
Apache License 2.0 (see `LICENSE`). Copyright 2026 Odinzen LLC.
