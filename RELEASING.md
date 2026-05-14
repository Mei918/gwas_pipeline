# Releasing `gwas-pipeline`

This repository uses a minimal tag-driven GitHub Release flow.

## Before tagging

1. Update the version in:
   - `pyproject.toml`
   - `src/gwas_pipeline/__init__.py`
2. Add an entry to `CHANGELOG.md`
3. Run the local checks you care about:

```bash
PYTHONPATH=src python3 -m gwas_pipeline --help
PYTHONPATH=src python3 -m gwas_pipeline doctor --help
```

4. Commit and push `main`

## Create a release tag

```bash
git tag v0.3.1
git push origin v0.3.1
```

## What happens next

Pushing a `v*` tag triggers `.github/workflows/release.yml`, which:

- installs the package with `.[dev]`
- builds `sdist` and `wheel`
- creates a GitHub Release
- uploads `dist/*` as release assets

## Notes

- The CI workflow remains separate from the release workflow.
- If you want PyPI publishing later, it can be added on top of this without
  changing the tag flow.
