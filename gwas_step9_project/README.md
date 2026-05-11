# gwas-step9-project

Standalone Python pipeline for GWAS step 9:

- sample-type geographic map visualization
- haplotype geographic map visualization
- haplotype proportion versus latitude regression summary

This project follows the workflow described in the blog's ninth step and
keeps the analysis separate from other projects.

## Input formats

### Sample-type map input

Headerless CSV with 4 columns:

- `Type`
- `Latitude`
- `Longitude`
- `City`

Example:

```text
Landrace,40.71048785,-74.00123991,New York
Wild,51.51669898,-0.130833151,London
```

### Haplotype map input

Headerless CSV with 6 columns:

- `Sample`
- `Haplotype`
- `Type`
- `Latitude`
- `Longitude`
- `City`

Example:

```text
S317,Hap3,Wild,-34.61284826,-58.37733518,Buenos Aires
S321,Hap2,Wild,-33.93340891,18.41515763,Cape Town
```

## Run

```bash
python3 gwas_step9_project/step9_geo_visualization.py \
  --samples /path/to/samples.csv \
  --hap-samples /path/to/hap_samples.csv \
  --output-dir /path/to/gwas_step9_output
```

## Optional regional zoom

```bash
python3 gwas_step9_project/step9_geo_visualization.py \
  --samples /path/to/samples.csv \
  --hap-samples /path/to/hap_samples.csv \
  --output-dir /path/to/gwas_step9_output \
  --xlim -10 50 \
  --ylim 20 70
```

## Main outputs

- `sample_type_counts_by_site.tsv`
- `haplotype_counts_by_site.tsv`
- `latitude_regression.tsv`
- `plots/sample_type_map.svg`
- `plots/haplotype_map.svg`
- `plots/latitude_regression.svg`
- `logs/pipeline.log`
- `state/pipeline_state.json`

## Notes

- Maps are rendered as lightweight SVGs with pie charts placed in geographic
  coordinates.
- Pie size is proportional to the square root of the number of samples at each
  site, similar to the blog's visual logic.
- The latitude regression uses simple linear regression computed in pure Python.
- Reruns skip completed steps unless `--force` is used.
