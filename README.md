# gwas-pipeline

`gwas-pipeline` is an installable Python package that wraps a 9-step GWAS
tutorial workflow behind a single CLI.

It keeps the original step-by-step tutorial layout, but also provides a
package-first interface:

```bash
gwas-pipeline step1 ...
gwas-pipeline step2 ...
...
gwas-pipeline step9 ...
```

## What is included

- packaged implementations for `step1` to `step9`
- thin legacy wrapper scripts in `gwas_step*_project/`
- server-focused runbooks and command examples
- a `doctor` command for checking environment dependencies

## Install

```bash
pip install -e .
```

## Quick Start

```bash
gwas-pipeline --help
gwas-pipeline doctor --profile plink_env
gwas-pipeline step1 --help
gwas-pipeline step8 --help
python -m gwas_pipeline step9 --help
```

## Repository Layout

- `src/gwas_pipeline/`
  Packaged step implementations and CLI entrypoints.
- `gwas_step1_project/` to `gwas_step9_project/`
  Legacy-compatible wrappers and per-step tutorial READMEs.
- `PACKAGE_COMMANDS.md`
  Copy-paste command examples for all steps.
- `SERVER_GWAS_RUNBOOK.md`
  Short command-first server runbook.
- `server_setup.md`
  Conda environment, install, and deployment notes.

## Tested Status

The tutorial workflow has been tested step-by-step on server.

- `step1`
  Real public FASTQ subset plus real `chr22.fa`
- `step2`
  Real chr22-region VCF converted to PLINK
- `step3`
  Real GFF plus synthetic GWAS result for logic validation
- `step4`
  Synthetic GWAS result for plotting validation
- `step5` to `step8`
  Real public human chr22 / COMT-region inputs
- `step9`
  Real sample geography tables

## Notes

- This repository is intentionally limited to the GWAS tutorial package.
- The separate `gwas_postgwas_tools` project is not included here.
- `gwas-pipeline` is the primary CLI; the old step scripts remain for
  compatibility with the original tutorial structure.
