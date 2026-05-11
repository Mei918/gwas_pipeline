# gwas-step2-project

Standalone Python pipeline for GWAS step 2:

- prepare phenotype files
- build GRM and PCA covariates
- run GCTA MLM GWAS
- generate Manhattan plot and QQ plot as SVG
- optionally write and run an rMVP R script

This project follows the workflow described in the blog's second step and
keeps the analysis separate from other projects.

## Inputs for GCTA mode

- PLINK prefix:
  - `cohort.bed`
  - `cohort.bim`
  - `cohort.fam`
- phenotype file like:
  - `S1<TAB>34.54`
  - `S10<TAB>35.95`

## Run GCTA MLM

```bash
python3 gwas_step2_project/step2_gwas_analysis.py \
  --mode gcta \
  --bfile /path/to/cohort \
  --phenotype /path/to/phenotype_of_NessMonster.txt \
  --output-dir /path/to/gwas_step2_output \
  --pca-components 5
```

## Run rMVP

```bash
python3 gwas_step2_project/step2_gwas_analysis.py \
  --mode rmvp \
  --vcf /path/to/filtersnp_fixed.vcf \
  --phenotype /path/to/phenotype_of_NessMonster.txt \
  --output-dir /path/to/gwas_step2_output
```

## Main outputs in GCTA mode

- `phenotype.txt`
- `gwas_pca.eigenvec`
- `kinship_matrix.grm.bin`
- `gwas_mlm_results.mlma`
- `plots/manhattan.svg`
- `plots/qq.svg`
- `logs/pipeline.log`
- `state/pipeline_state.json`

## Notes

- The script updates the phenotype-linked FAM as a copied file in the output
  directory and does not overwrite the original input `cohort.fam`.
- Reruns will skip finished steps unless `--force` is used.
