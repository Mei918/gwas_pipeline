# Server Setup

This document covers server environment creation, dependency installation,
package installation, and deployment flow for the packaged GWAS tutorial
pipeline.

Server login:

```bash
ssh -p 13109 y413109@ssh.sxqtx.com
```

## 1. Conda Environments

Suggested environments:
- `plink_env`: alignment, variant, PLINK, and GCTA toolchain
- `omicverse`: plotting and Python-heavy downstream steps

List existing environments:

```bash
conda env list
```

Create if needed:

```bash
conda create -n plink_env python=3.10 -y
conda create -n omicverse python=3.10 -y
```

## 2. Install CLI Dependencies

### plink_env

```bash
conda activate plink_env
conda install -c bioconda plink bcftools vcftools bwa samtools fastp gatk4 gcta -y
```

Verify:

```bash
which plink
which bcftools
which vcftools
which bwa
which samtools
which fastp
which gatk
which gcta64
```

### omicverse

```bash
conda activate omicverse
conda install pandas numpy matplotlib networkx pyyaml pip -y
```

Verify:

```bash
python3 - <<'PY'
mods = ["pandas", "numpy", "matplotlib", "networkx", "yaml"]
for m in mods:
    try:
        __import__(m)
        print(m, "OK")
    except Exception as e:
        print(m, "FAIL", e)
PY
```

## 3. Deploy The Repository

Copy code from local machine:

```bash
scp -P 13109 -r "/Users/may/Documents/New project" y413109@ssh.sxqtx.com:/home/y413109/gwas_project_src
```

Or sync updates:

```bash
rsync -av -e "ssh -p 13109" "/Users/may/Documents/New project/" y413109@ssh.sxqtx.com:/home/y413109/gwas_project_src/
```

If you already work directly inside `~/gwas_project`, keep the repository there instead.

## 4. Install The Package

From the repository root on the server:

```bash
cd ~/gwas_project_src
conda activate omicverse
pip install -e .
```

If you prefer installing from `plink_env`, that also works:

```bash
cd ~/gwas_project_src
conda activate plink_env
pip install -e .
```

Verify package entrypoints:

```bash
gwas-pipeline --help
python -m gwas_pipeline --help
```

Verify a representative step:

```bash
gwas-pipeline step3 --help
gwas-pipeline step8 --help
```

## 5. Runtime Guidance

Recommended environment by step:
- Step 1: `plink_env`
- Step 2: `plink_env`
- Step 3: `omicverse`
- Step 4: `omicverse`
- Step 5: `omicverse`
- Step 6: `omicverse`
- Step 7: `omicverse`
- Step 8: `plink_env`
- Step 9: `omicverse`

If the shell cannot find `gwas-pipeline` after install:

```bash
python -m gwas_pipeline --help
```

## 6. Common Checks

Show package metadata:

```bash
python - <<'PY'
import gwas_pipeline
print(gwas_pipeline.__version__)
PY
```

Check old wrappers still work:

```bash
python3 gwas_step1_project/step1_fastq_to_vcf.py --help
python3 gwas_step8_project/step8_pi_tajima.py --help
```

## 7. Main Workflow Documents

- [SERVER_GWAS_RUNBOOK.md](./SERVER_GWAS_RUNBOOK.md)
- [PACKAGE_COMMANDS.md](./PACKAGE_COMMANDS.md)
