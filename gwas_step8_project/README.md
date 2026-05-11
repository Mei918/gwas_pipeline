# gwas-step8-project

Standalone Python pipeline for GWAS step 8:

- extract a target region from a VCF
- split samples by population from a classification CSV
- run `vcftools` to compute nucleotide diversity (pi) and Tajima's D
- merge per-group outputs
- draw a two-panel SVG figure for pi and Tajima's D

This project follows the workflow described in the blog's eighth step and
keeps the analysis separate from other projects.

## Inputs

- whole-genome or regional VCF
- classification CSV such as:
  - `Sample,Group`
  - `S1,Cultivar`
  - `S10,Landrace`

## Run

```bash
python3 gwas_step8_project/step8_pi_tajima.py \
  --vcf /path/to/combine_gatk.vcf.gz \
  --region Chr3:66400000-66411250 \
  --class-file /path/to/clss.csv \
  --output-dir /path/to/gwas_step8_output \
  --window-size 100 \
  --step-size 100
```

## Main outputs

- `gene_region.vcf.gz`
- `groups/*.samples.txt`
- `groups/*.vcf.gz`
- `stats/*_pi.windowed.pi`
- `stats/*_tajima.Tajima.D`
- `merged/pi_merged.tsv`
- `merged/tajima_merged.tsv`
- `plots/pi_tajima.svg`
- `logs/pipeline.log`
- `state/pipeline_state.json`

## Notes

- Significant Tajima's D blocks are highlighted when `TajimaD < -2`.
- By default the grouping column is the second column in the class file, but you
  can override it with `--group-column`.
- Reruns skip completed steps unless `--force` is used.
