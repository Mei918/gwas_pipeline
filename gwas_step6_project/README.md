# gwas-step6-project

Standalone Python pipeline for GWAS step 6:

- assign haplotypes from a regional VCF
- label haplotypes as Hap1..HapN by frequency
- merge haplotypes with sample metadata
- summarize haplotype composition by metadata group
- draw a proportional stacked SVG plot

This project follows the workflow described in the blog's sixth step and
keeps the analysis separate from other projects.

## Input files

- regional VCF, for example `gene.vcf`
- sample metadata CSV, for example:
  - `Sample,Region,Breeding_time`

## Run

```bash
python3 gwas_step6_project/step6_haplotype_analysis.py \
  --vcf /path/to/gene.vcf \
  --resource /path/to/resource.csv \
  --output-dir /path/to/gwas_step6_output \
  --group-column Region
```

## Main outputs

- `gene.csv`
- `haplotype_sample_table.tsv`
- `haplotype_group_summary.tsv`
- `plots/haplotype_by_group.svg`
- `logs/pipeline.log`
- `state/pipeline_state.json`

## Notes

- Haplotypes are ranked by sample count, so the most frequent haplotype is
  always `Hap1`.
- Missing genotypes such as `./.` are kept in the haplotype string.
- `--group-column` can be `Region`, `Breeding_time`, or any other metadata
  column present in the resource file.
- Reruns skip completed steps unless `--force` is used.
