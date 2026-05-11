# gwas-step3-project

Standalone Python pipeline for GWAS step 3:

- extract significant loci from GWAS results
- convert GFF gene records into gene-plus-promoter BED intervals
- annotate significant SNPs to overlapping genes

This project follows the workflow described in the blog's third step and
keeps the analysis separate from other projects.

## Supported inputs

- GWAS result table from tools such as PLINK or GCTA
- GFF file containing `gene` features

## Run

```bash
python3 gwas_step3_project/step3_candidate_gene_extraction.py \
  --gwas /path/to/fallgwas_glm.ADD.txt \
  --gff /path/to/nessmonster.gff \
  --output-dir /path/to/gwas_step3_output \
  --p-threshold 1e-5 \
  --chr-column CHR \
  --pos-column BP
```

## Main outputs

- `significant_hits.tsv`
- `genes_with_promoter.bed`
- `significant_hits.with_gene.tsv`
- `candidate_gene_summary.tsv`
- `logs/pipeline.log`
- `state/pipeline_state.json`

## Notes

- By default, promoters extend 2000 bp upstream for `+` strand genes and
  2000 bp downstream for `-` strand genes.
- The script can optionally prefix numeric chromosome values, for example
  converting `1` to `Chr1` to match BED naming.
- Reruns skip completed steps unless `--force` is used.
