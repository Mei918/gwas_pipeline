# Release Checklist

## Code

- [ ] `python -m py_compile src/gwas_pipeline/*.py src/gwas_pipeline/steps/*.py`
- [ ] `PYTHONPATH=src python3 -m gwas_pipeline --help`
- [ ] `PYTHONPATH=src python3 -m gwas_pipeline doctor --help`
- [ ] `PYTHONPATH=src python3 -m gwas_pipeline step1 --help`
- [ ] `PYTHONPATH=src python3 -m gwas_pipeline step8 --help`

## Packaging

- [ ] `pip install -e .`
- [ ] `gwas-pipeline --help`
- [ ] `gwas-pipeline doctor --profile plink_env`
- [ ] `gwas-pipeline doctor --profile omicverse`

## Documentation

- [ ] README describes only the GWAS tutorial package
- [ ] `PACKAGE_COMMANDS.md` matches current CLI flags
- [ ] `SERVER_GWAS_RUNBOOK.md` uses `gwas-pipeline stepX` commands
- [ ] `server_setup.md` reflects current conda environment setup
- [ ] `CHANGELOG.md` updated for this release

## Server Validation

- [ ] `step1` tested with real FASTQ subset and real reference
- [ ] `step2` tested on real chr22-region PLINK input
- [ ] `step5` to `step8` tested on real public human chr22 inputs
- [ ] `step9` tested on real sample geography tables

## GitHub

- [ ] clean `git status`
- [ ] commit created
- [ ] branch pushed
- [ ] repository homepage renders README clearly
