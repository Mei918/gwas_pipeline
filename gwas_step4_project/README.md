# gwas-step4-project

Standalone Python pipeline for GWAS step 4:

- extract a local genomic region from GWAS results
- draw a local Manhattan plot in SVG
- optionally overlay gene models from a BED file
- optionally run `LDBlockShow` for LD visualization

This project follows the workflow described in the blog's fourth step and
keeps the analysis separate from other projects.

## Inputs

- GWAS result table
- optional gene BED file
- optional VCF file for LD analysis

## Run local Manhattan only

```bash
python3 gwas_step4_project/step4_local_manhattan_ld.py \
  --gwas /path/to/gwas_mlm_results.txt \
  --output-dir /path/to/gwas_step4_output \
  --chr 3 \
  --start 63837500 \
  --end 63900000 \
  --chr-column CHR \
  --pos-column BP \
  --p-column P
```

## Run local Manhattan with gene track and LD

```bash
python3 gwas_step4_project/step4_local_manhattan_ld.py \
  --gwas /path/to/gwas_mlm_results.txt \
  --gene-bed /path/to/genes_only.bed \
  --vcf /path/to/filtered_snp.vcf \
  --ldblockshow /path/to/LDBlockShow \
  --output-dir /path/to/gwas_step4_output \
  --chr 8 \
  --start 23520000 \
  --end 23530000 \
  --add-chr-prefix \
  --chr-prefix Chr
```

## Main outputs

- `region_gwas.tsv`
- `plots/local_manhattan.svg`
- `plots/local_manhattan_with_genes.svg`
- `ld/test_output.*` from `LDBlockShow` when enabled
- `logs/pipeline.log`
- `state/pipeline_state.json`

## Notes

- The script accepts either whitespace-delimited or tab-delimited GWAS tables.
- Gene BED is expected as:
  - `CHR  start0  end0  gene_id  score  strand`
- Numeric chromosome values can be normalized with `--add-chr-prefix`.
- Reruns skip completed steps unless `--force` is used.
