# Changelog

## 0.3.1 - 2026-05-14

- Added a PyPI publishing workflow for manual GitHub Actions releases.
- Tightened CI to install `.[dev]` and build distributions during validation.
- Added `twine` to development dependencies.
- Documented PyPI installation in `README.md`.
- Expanded `RELEASING.md` with local build, `twine check`, TestPyPI, and PyPI
  upload commands.

## 0.3.0 - 2026-05-11

- Created a standalone `gwas_pipeline` package for the GWAS tutorial workflow.
- Added packaged CLI entrypoints for `step1` through `step9`.
- Added `gwas-pipeline doctor` for environment dependency checks.
- Preserved legacy `gwas_step*_project` scripts as thin compatibility wrappers.
- Added server-oriented documentation:
  - `PACKAGE_COMMANDS.md`
  - `SERVER_GWAS_RUNBOOK.md`
  - `server_setup.md`
- Folded server-tested compatibility fixes into packaged modules:
  - `step1`: read-group repair, force cleanup for reference indexes, graceful
    skip when filtered VCF has no variants for PLINK conversion
  - `step8`: VCFv4.2 normalization for `vcftools`, no `--TajimaD-step`, and
    fallback for missing `BIN_END`
