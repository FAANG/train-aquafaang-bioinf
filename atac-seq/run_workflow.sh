#!/bin/bash

# bash strict mode, i.e. stop script on first error
set -euo pipefail

# need snakemake
module load snakemake

# need samtools
module load SAMtools


snakemake \
  --cores 16 \
  --printshellcmds \
  -s workflow/Snakefile \
  --configfile config/workflow_config.yaml \
  --keep-going

# Dry-run:
# snakemake -n -j1 --printshellcmds -s workflow/Snakefile --configfile config/workflow_config.yaml
