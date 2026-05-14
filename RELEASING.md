# Releasing `gwas-pipeline`

This repository uses a minimal tag-driven GitHub Release flow and supports
manual PyPI publishing with `twine`.

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

## Build and verify distributions locally

```bash
python -m pip install -e ".[dev]"
python -m build
python -m twine check dist/*
```

## Upload to PyPI manually

```bash
python -m twine upload dist/*
```

When prompted by `twine`, use:

- username: `__token__`
- password: your PyPI API token

If you want to test the publish flow first:

```bash
python -m twine upload --repository testpypi dist/*
```

## Create a release tag

```bash
git tag v0.3.1
git push origin v0.3.1
```

## Copy-paste release sequence

```bash
cd "/Users/may/Documents/New project/gwas_tutorial_repo"
python -m pip install -e ".[dev]"
python -m build
python -m twine check dist/*
python -m twine upload dist/*
git add pyproject.toml src/gwas_pipeline/__init__.py CHANGELOG.md README.md RELEASING.md
git commit -m "Release 0.3.1"
git push origin main
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
