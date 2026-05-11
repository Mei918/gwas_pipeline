# gwas-pipeline

Installable Python package for the GWAS tutorial workflow, with packaged
`step1` through `step9` modules and matching legacy wrapper scripts.

## Install

```bash
pip install -e .
```

## CLI

```bash
gwas-pipeline --help
gwas-pipeline doctor --profile plink_env
gwas-pipeline step1 --help
gwas-pipeline step8 --help
python -m gwas_pipeline step9 --help
```

## Repository Layout

- `src/gwas_pipeline/`: packaged implementation of steps `1` to `9`
- `gwas_step*_project/`: thin compatibility wrappers plus per-step README files
- `PACKAGE_COMMANDS.md`: package command examples
- `SERVER_GWAS_RUNBOOK.md`: server runbook focused on package commands
- `server_setup.md`: conda environment and deployment notes

## Tested Workflow Status

The tutorial steps were validated on server with a mix of real public human
chr22 test data and minimal synthetic inputs:

- `step1`: real FASTQ subset + real `chr22.fa`
- `step2`: real VCF-derived PLINK test set
- `step3`: real GFF plus synthetic GWAS result for logic validation
- `step4`: synthetic GWAS result for plotting validation
- `step5` to `step8`: real public human chr22/COMT-region inputs
- `step9`: real sample geography tables

## Notes

- The separate `gwas_postgwas_tools` project is intentionally not included in
  this repository.
- The packaged CLI is `gwas-pipeline`; legacy wrapper scripts remain for
  compatibility with the original tutorial layout.
