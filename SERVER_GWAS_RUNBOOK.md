# Server GWAS Runbook

```bash
ssh -p 13109 y413109@ssh.sxqtx.com
cd ~/gwas_project
pip install -e .
gwas-pipeline doctor
```

```bash
mkdir -p ~/gwas_project/human_chr22_test
cd ~/gwas_project/human_chr22_test
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi
wget https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
```

```bash
mkdir -p ~/gwas_project/reference_chr22
cd ~/gwas_project/reference_chr22
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz
gunzip chr22.fa.gz
```

```bash
mkdir -p ~/gwas_project/human_fastq_test
cd ~/gwas_project/human_fastq_test
wget -O ena_file_report.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR622461&result=read_run&fields=fastq_ftp"
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR622/SRR622461/SRR622461_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR622/SRR622461/SRR622461_2.fastq.gz
gzip -t SRR622461_1.fastq.gz
gzip -t SRR622461_2.fastq.gz
gunzip -c SRR622461_1.fastq.gz | awk 'NR<=400000' | gzip > test_R1.fastq.gz
gunzip -c SRR622461_2.fastq.gz | awk 'NR<=400000' | gzip > test_R2.fastq.gz
cp test_R1.fastq.gz Z1_UDI601.R1.fastq.gz
cp test_R2.fastq.gz Z1_UDI601.R2.fastq.gz
cp test_R1.fastq.gz Z2_UDI601.R1.fastq.gz
cp test_R2.fastq.gz Z2_UDI601.R2.fastq.gz
```

```bash
cd ~/gwas_project/human_chr22_test
awk 'BEGIN{OFS="\t"} NR>1 && $3=="AFR" && c1<4 {print $1,"AFR"; c1++}
     NR>1 && $3=="EUR" && c2<4 {print $1,"EUR"; c2++}
     NR>1 && $3=="EAS" && c3<4 {print $1,"EAS"; c3++}' \
integrated_call_samples_v3.20130502.ALL.panel > selected_samples.tsv
cut -f1 selected_samples.tsv > selected_samples.txt
awk 'BEGIN{print "Sample,Group"} {print $1","$2}' selected_samples.tsv > class.csv
awk 'BEGIN{print "Sample,Region,Breeding_time"} {print $1","$2",TestEra"}' selected_samples.tsv > resource.csv
bcftools view ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -r 22:19941371-19969975 -Oz -o COMT_region_all.vcf.gz
bcftools index -t COMT_region_all.vcf.gz
bcftools view -S selected_samples.txt COMT_region_all.vcf.gz -Oz -o COMT_region_12samples.vcf.gz
bcftools index -t COMT_region_12samples.vcf.gz
bcftools view COMT_region_12samples.vcf.gz -Ov -o COMT_region_12samples.vcf
gunzip -c Homo_sapiens.GRCh38.115.gff3.gz > Homo_sapiens.GRCh38.115.gff3
gunzip -c Homo_sapiens.GRCh38.cds.all.fa.gz > Homo_sapiens.GRCh38.cds.all.fa
cat > mini_gwas.txt <<'EOF'
CHR	SNP	BP	A1	TEST	NMISS	BETA	STAT	P
22	rs1	19941400	A	ADD	12	0.25	4.50	1e-8
22	rs2	19941520	C	ADD	12	-0.32	-4.10	2e-7
22	rs3	19941680	G	ADD	12	0.41	5.20	6e-6
22	rs4	19941810	T	ADD	12	-0.29	-4.80	8e-6
22	rs5	19942000	A	ADD	12	0.05	0.80	0.42
22	rs6	19942200	C	ADD	12	-0.03	-0.55	0.58
EOF
plink --vcf COMT_region_12samples.vcf --make-bed --out cohort
awk 'BEGIN{OFS="\t"} {print $2, 10+NR}' cohort.fam > phenotype_of_test.txt
```

```bash
conda activate plink_env
gwas-pipeline step1 --fastq-dir ~/gwas_project/human_fastq_test --reference ~/gwas_project/reference_chr22/chr22.fa --output-dir ~/gwas_project/human_fastq_test/output_real_chr22 --fastp-threads 2 --bwa-threads 2 --sort-threads 2 --maf 0 --max-missing 1.0 --force
gwas-pipeline step2 --mode gcta --bfile ~/gwas_project/human_chr22_test/cohort --phenotype ~/gwas_project/human_chr22_test/phenotype_of_test.txt --output-dir ~/gwas_project/human_chr22_test/test_step2_output --pca-components 3
gwas-pipeline step8 --vcf ~/gwas_project/human_chr22_test/COMT_region_all.vcf.gz --region 22:19941371-19969975 --class-file ~/gwas_project/human_chr22_test/class.csv --output-dir ~/gwas_project/human_chr22_test/test_step8_output_fix2 --sample-column Sample --group-column Group --window-size 100 --step-size 100
```

```bash
conda activate omicverse
gwas-pipeline step3 --gwas ~/gwas_project/human_chr22_test/mini_gwas.txt --gff ~/gwas_project/human_chr22_test/Homo_sapiens.GRCh38.115.gff3 --output-dir ~/gwas_project/human_chr22_test/test_step3_output_fix --p-threshold 1e-5 --p-column P --chr-column CHR --pos-column BP
gwas-pipeline step4 --gwas ~/gwas_project/human_chr22_test/mini_gwas.txt --output-dir ~/gwas_project/human_chr22_test/test_step4_output --chr 22 --start 19941371 --end 19969975 --chr-column CHR --pos-column BP --p-column P
gwas-pipeline step5 --region-vcf ~/gwas_project/human_chr22_test/COMT_region_12samples.vcf --gff ~/gwas_project/human_chr22_test/Homo_sapiens.GRCh38.115.gff3 --cds-fasta ~/gwas_project/human_chr22_test/Homo_sapiens.GRCh38.cds.all.fa --output-dir ~/gwas_project/human_chr22_test/test_step5_output
gwas-pipeline step6 --vcf ~/gwas_project/human_chr22_test/COMT_region_12samples.vcf --resource ~/gwas_project/human_chr22_test/resource.csv --output-dir ~/gwas_project/human_chr22_test/test_step6_output --group-column Region --group-order AFR,EUR,EAS
gwas-pipeline step7 --vcf ~/gwas_project/human_chr22_test/COMT_region_12samples.vcf --resource ~/gwas_project/human_chr22_test/resource.csv --output-dir ~/gwas_project/human_chr22_test/test_step7_output --group-column Region --group-order AFR,EUR,EAS --network-k 2.0 --network-seed 263
gwas-pipeline step9 --samples ~/gwas_project/samples.csv --hap-samples ~/gwas_project/hap_samples.csv --output-dir ~/gwas_project/test_step9_output
```
