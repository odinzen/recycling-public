#!/usr/bin/env bash
# Reproduce every calculated number, table and figure of the titanium transus paper.
# Needs the reduced Ti-O and Ti-N databases in data/ (available from the corresponding author) and pycalphad's
# bundled Ghosh 2002 Cr-Ti-V database.
set -e
python ti_validation.py            # pure Ti and Ti-V transus      -> ti_validation_results.json, Fig. 1
python ti_oeq_coefficient_test.py  # O and N transus, N:O weight   -> ti_oeq_coefficient_results.json, Fig. 2
python ti_partition.py             # tie lines, partition coeff.   -> ti_partition_results.json, Fig. 3
python ti_transfer_entropy.py      # enthalpy/entropy of transfer  -> ti_transfer_entropy_results.json, Fig. 4
