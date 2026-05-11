# gwas-step5-project

Standalone Python pipeline for GWAS step 5:

- optionally extract a target region from a VCF with `bcftools`
- convert regional VCF into a sample-by-variant CSV table
- parse GFF and CDS FASTA
- detect whether coding variants are synonymous or nonsynonymous
- rewrite variant columns with amino-acid effect annotations

This project follows the workflow described in the blog's fifth step and
keeps the analysis separate from other projects.

## Run from an existing regional VCF

```bash
python3 gwas_step5_project/step5_amino_acid_mutation.py \
  --region-vcf /path/to/gene_region.vcf \
  --gff /path/to/genome.gff \
  --cds-fasta /path/to/cds.fasta \
  --output-dir /path/to/gwas_step5_output
```

## Run from a whole VCF plus region

```bash
python3 gwas_step5_project/step5_amino_acid_mutation.py \
  --vcf /path/to/combine_gatk.vcf.gz \
  --region Chr4:49881219-49887065 \
  --gff /path/to/genome.gff \
  --cds-fasta /path/to/cds.fasta \
  --output-dir /path/to/gwas_step5_output
```

## Main outputs

- `gene_region.vcf`
- `gene_region.csv`
- `gene_region_mut.csv`
- `variant_effect_summary.tsv`
- `logs/pipeline.log`
- `state/pipeline_state.json`

## Notes

- Variants that overlap CDS and change amino acids are labeled like
  `NonSynonymous_Q->K`.
- Variants overlapping CDS but not changing amino acid are labeled
  `Synonymous`.
- Variants without a CDS effect keep their original column name.
- Reruns skip finished steps unless `--force` is used.
