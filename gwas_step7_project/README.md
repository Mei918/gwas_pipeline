# gwas-step7-project

Standalone Python pipeline for GWAS step 7:

- parse a regional VCF into haplotype sequences
- collapse samples into unique haplotypes
- compute pairwise haplotype distances
- build a minimum spanning haplotype network
- draw a pie-chart haplotype network as SVG

This project follows the workflow described in the blog's seventh step and
keeps the analysis separate from other projects.

## Inputs

- regional VCF file
- resource CSV with at least:
  - `Sample`
  - a grouping column such as `Region`

## Run

```bash
python3 gwas_step7_project/step7_haplotype_network.py \
  --vcf /path/to/gene.vcf \
  --resource /path/to/resource.csv \
  --output-dir /path/to/gwas_step7_output \
  --group-column Region \
  --network-k 3.2 \
  --network-seed 263
```

## Main outputs

- `haplotypes.tsv`
- `haplotype_samples.tsv`
- `haplotype_distances.tsv`
- `mst_edges.tsv`
- `hap_network.svg`
- `logs/pipeline.log`
- `state/pipeline_state.json`

## Notes

- The most frequent haplotype is labeled `Hap1`.
- The network layout is controlled by `--network-k`, `--network-seed`, and
  `--network-iterations`.
- The plotting step uses `networkx`, `matplotlib`, `numpy`, and `pandas`.
- Reruns skip completed steps unless `--force` is used.
