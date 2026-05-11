# Package Commands

Install editable package from the repository root:

```bash
pip install -e .
```

Primary entrypoints after install:

```bash
gwas-pipeline --help
python -m gwas_pipeline --help
```

## GWAS Tutorial Steps

### Step 1

```bash
gwas-pipeline step1 \
  --fastq-dir ~/gwas_project/human_fastq_test \
  --reference ~/gwas_project/reference_chr22/chr22.fa \
  --output-dir ~/gwas_project/human_fastq_test/output_real_chr22 \
  --fastp-threads 2 \
  --bwa-threads 2 \
  --sort-threads 2 \
  --maf 0 \
  --max-missing 1.0 \
  --force
```

### Step 2

```bash
gwas-pipeline step2 \
  --mode gcta \
  --bfile ~/gwas_project/human_chr22_test/cohort \
  --phenotype ~/gwas_project/human_chr22_test/phenotype_of_test.txt \
  --output-dir ~/gwas_project/human_chr22_test/test_step2_output \
  --pca-components 3
```

### Step 3

```bash
gwas-pipeline step3 \
  --gwas ~/gwas_project/human_chr22_test/mini_gwas.txt \
  --gff ~/gwas_project/human_chr22_test/Homo_sapiens.GRCh38.115.gff3 \
  --output-dir ~/gwas_project/human_chr22_test/test_step3_output_fix \
  --p-threshold 1e-5 \
  --p-column P \
  --chr-column CHR \
  --pos-column BP
```

### Step 4

```bash
gwas-pipeline step4 \
  --gwas ~/gwas_project/human_chr22_test/mini_gwas.txt \
  --output-dir ~/gwas_project/human_chr22_test/test_step4_output \
  --chr 22 \
  --start 19941371 \
  --end 19969975 \
  --chr-column CHR \
  --pos-column BP \
  --p-column P
```

### Step 5

```bash
gwas-pipeline step5 \
  --region-vcf ~/gwas_project/human_chr22_test/COMT_region_12samples.vcf \
  --gff ~/gwas_project/human_chr22_test/Homo_sapiens.GRCh38.115.gff3 \
  --cds-fasta ~/gwas_project/human_chr22_test/Homo_sapiens.GRCh38.cds.all.fa \
  --output-dir ~/gwas_project/human_chr22_test/test_step5_output
```

### Step 6

```bash
gwas-pipeline step6 \
  --vcf ~/gwas_project/human_chr22_test/COMT_region_12samples.vcf \
  --resource ~/gwas_project/human_chr22_test/resource.csv \
  --output-dir ~/gwas_project/human_chr22_test/test_step6_output \
  --group-column Region \
  --group-order AFR,EUR,EAS
```

### Step 7

```bash
gwas-pipeline step7 \
  --vcf ~/gwas_project/human_chr22_test/COMT_region_12samples.vcf \
  --resource ~/gwas_project/human_chr22_test/resource.csv \
  --output-dir ~/gwas_project/human_chr22_test/test_step7_output \
  --group-column Region \
  --group-order AFR,EUR,EAS \
  --network-k 2.0 \
  --network-seed 263
```

### Step 8

```bash
gwas-pipeline step8 \
  --vcf ~/gwas_project/human_chr22_test/COMT_region_all.vcf.gz \
  --region 22:19941371-19969975 \
  --class-file ~/gwas_project/human_chr22_test/class.csv \
  --output-dir ~/gwas_project/human_chr22_test/test_step8_output_fix2 \
  --sample-column Sample \
  --group-column Group \
  --window-size 100 \
  --step-size 100
```

### Step 9

```bash
gwas-pipeline step9 \
  --samples ~/gwas_project/samples.csv \
  --hap-samples ~/gwas_project/hap_samples.csv \
  --output-dir ~/gwas_project/test_step9_output
```

## Legacy Compatibility

The old step scripts still work and now forward to the packaged modules:

```bash
python3 gwas_step8_project/step8_pi_tajima.py --help
python3 gwas_step1_project/step1_fastq_to_vcf.py --help
```
