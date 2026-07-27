#!/usr/bin/env bash
# Reproduce the computed results and figures for the recycled aerospace Ti study.
# Runs against the reduced metallic TDBs in data/ and pycalphad's bundled Ghosh 2002.
set -e
python ti_validation.py            # pure-Ti + Ti-V transus   -> ti_validation_results.json, figure5
python ti_o_transus.py             # Ti-O transus shift       -> ti_o_transus_results.json, figure6
python ti_o_partition.py           # Ti-O partition, 2-phase  -> ti_o_partition_results.json, figure7
python ti_oeq_coefficient_test.py  # N vs O potency           -> ti_oeq_coefficient_results.json, figure8
python ti_transus_figures.py       # schematics               -> figure1..4
echo "Done. See the *_results.json files, figures/, and CIRCULARITY_AUDIT.md."
