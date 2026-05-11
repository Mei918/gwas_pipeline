# gwas-step1-project

Standalone Python pipeline for GWAS step 1:

- scan paired-end `fastq.gz` files automatically
- run QC with `fastp`
- build reference indexes
- align reads with `bwa`
- sort and index BAM with `samtools`
- call per-sample GVCF with `gatk HaplotypeCaller`
- combine GVCFs and genotype cohort VCF
- filter VCF with `vcftools`
- convert VCF to PLINK with `plink`

This version adds:

- automatic sample discovery
- log files
- resume support
- skip completed steps

## Expected FASTQ naming

The scanner expects filenames like:

- `Z1_UDI601.R1.fastq.gz`
- `Z1_UDI601.R2.fastq.gz`

It groups files by the shared prefix before `.R1.fastq.gz` or `.R2.fastq.gz`.

## Run

```bash
python3 gwas_step1_project/step1_fastq_to_vcf.py \
  --fastq-dir /path/fastq \
  --reference /path/zm1.fasta \
  --output-dir /path/step1_output \
  --fastp-threads 8 \
  --bwa-threads 32 \
  --sort-threads 8 \
  --maf 0.05 \
  --max-missing 0.8
```

## Key outputs

- `logs/pipeline.log`
- `state/pipeline_state.json`
- `clean/`
- `qc/`
- `bam/`
- `gvcf/`
- `vcf/cohort.raw.vcf.gz`
- `vcf/cohort.gwas.recode.vcf`
- `plink/cohort.bed`

## Resume behavior

On rerun, the pipeline:

- reloads `state/pipeline_state.json`
- skips a step if its expected output already exists
- records completed steps by sample and by cohort stage

This makes it safe to restart after interruption.
