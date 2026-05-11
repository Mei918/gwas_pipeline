from __future__ import annotations

import argparse
import csv
import json
import logging
from collections import Counter, defaultdict
from pathlib import Path
import sys
from typing import Any


IUPAC = {
    ("A", "A"): "A",
    ("T", "T"): "T",
    ("C", "C"): "C",
    ("G", "G"): "G",
    ("A", "G"): "R",
    ("G", "A"): "R",
    ("C", "T"): "Y",
    ("T", "C"): "Y",
    ("G", "T"): "K",
    ("T", "G"): "K",
    ("A", "C"): "M",
    ("C", "A"): "M",
    ("C", "G"): "S",
    ("G", "C"): "S",
    ("A", "T"): "W",
    ("T", "A"): "W",
    (".", "."): "N",
}


class Step7Pipeline:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.output_dir = Path(args.output_dir).resolve()
        self.dirs = self._ensure_dirs()
        self.logger = self._build_logger()
        self.state_path = self.dirs["state"] / "pipeline_state.json"
        self.state = self._load_state()

    def run(self) -> None:
        self.logger.info("Starting GWAS step7 haplotype network workflow")
        haplotypes_path = self.extract_haplotypes()
        samples_path = self.expand_haplotype_samples(haplotypes_path)
        distances_path = self.compute_haplotype_distances(haplotypes_path)
        mst_path = self.compute_minimum_spanning_tree(distances_path)
        self.draw_network(haplotypes_path, samples_path, mst_path)
        self.logger.info("Step7 workflow finished successfully")

    def extract_haplotypes(self) -> Path:
        out_path = self.output_dir / "haplotypes.tsv"
        if self._should_skip_global("extract_haplotypes", [out_path]):
            self.logger.info("Skipping haplotype extraction")
            return out_path

        vcf_path = Path(self.args.vcf).resolve()
        if not vcf_path.exists():
            raise FileNotFoundError(f"Missing VCF file: {vcf_path}")

        samples: list[str] = []
        sample_seqs: dict[str, list[str]] = {}
        with vcf_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    samples = line.strip().split("\t")[9:]
                    sample_seqs = {sample: [] for sample in samples}
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 10:
                    continue
                alt_alleles = [cols[3]] + cols[4].split(",")
                for idx, gt_field in enumerate(cols[9:]):
                    gt = gt_field.split(":", 1)[0]
                    sample_seqs[samples[idx]].append(self._get_iupac(gt, alt_alleles))

        sample_to_seq = {sample: "".join(chars) for sample, chars in sample_seqs.items()}
        hap_to_samples: dict[str, list[str]] = defaultdict(list)
        for sample, seq in sample_to_seq.items():
            hap_to_samples[seq].append(sample)

        unique_haps = sorted(hap_to_samples.keys(), key=lambda seq: len(hap_to_samples[seq]), reverse=True)
        total_samples = len(samples)
        rows = []
        for idx, seq in enumerate(unique_haps, start=1):
            hap_name = f"Hap{idx}"
            count = len(hap_to_samples[seq])
            rows.append(
                {
                    "Haplotype": hap_name,
                    "Sequence": seq,
                    "Count": str(count),
                    "Frequency": f"{(count / total_samples) * 100:.2f}",
                    "Sample_List": ",".join(hap_to_samples[seq]),
                }
            )

        with out_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["Haplotype", "Sequence", "Count", "Frequency", "Sample_List"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rows)

        self.logger.info("Collapsed %d samples into %d haplotypes", total_samples, len(rows))
        self._mark_global_done("extract_haplotypes")
        return out_path

    def expand_haplotype_samples(self, haplotypes_path: Path) -> Path:
        out_path = self.output_dir / "haplotype_samples.tsv"
        if self._should_skip_global("expand_haplotype_samples", [out_path]):
            self.logger.info("Skipping haplotype sample expansion")
            return out_path

        resource_path = Path(self.args.resource).resolve()
        resource_rows = self._read_resource_rows(resource_path)
        resource_by_sample = {row["Sample"]: row for row in resource_rows}

        rows = []
        with haplotypes_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                hap = row["Haplotype"]
                seq = row["Sequence"]
                for sample in [item.strip() for item in row["Sample_List"].split(",") if item.strip()]:
                    resource = resource_by_sample.get(sample, {})
                    out_row = {"Sample": sample, "Haplotype": hap, "Sequence": seq}
                    out_row.update(resource)
                    rows.append(out_row)

        fieldnames = self._ordered_fieldnames(rows, ["Sample", "Haplotype", "Sequence"])
        with out_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

        self.logger.info("Expanded haplotypes into %d sample records", len(rows))
        self._mark_global_done("expand_haplotype_samples")
        return out_path

    def compute_haplotype_distances(self, haplotypes_path: Path) -> Path:
        out_path = self.output_dir / "haplotype_distances.tsv"
        if self._should_skip_global("compute_haplotype_distances", [out_path]):
            self.logger.info("Skipping haplotype distance table")
            return out_path

        hap_rows = self._read_tsv_rows(haplotypes_path)
        rows = []
        for i in range(len(hap_rows)):
            for j in range(i + 1, len(hap_rows)):
                seq1 = hap_rows[i]["Sequence"]
                seq2 = hap_rows[j]["Sequence"]
                distance = self._hamming_distance(seq1, seq2)
                rows.append(
                    {
                        "Haplotype1": hap_rows[i]["Haplotype"],
                        "Haplotype2": hap_rows[j]["Haplotype"],
                        "Distance": str(distance),
                    }
                )

        with out_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["Haplotype1", "Haplotype2", "Distance"], delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

        self.logger.info("Computed %d pairwise haplotype distances", len(rows))
        self._mark_global_done("compute_haplotype_distances")
        return out_path

    def compute_minimum_spanning_tree(self, distances_path: Path) -> Path:
        out_path = self.output_dir / "mst_edges.tsv"
        if self._should_skip_global("compute_minimum_spanning_tree", [out_path]):
            self.logger.info("Skipping minimum spanning tree")
            return out_path

        try:
            import networkx as nx
        except ImportError as exc:
            raise RuntimeError("networkx is required for step7 MST generation") from exc

        hap_rows = self._read_tsv_rows(self.output_dir / "haplotypes.tsv")
        distances = self._read_tsv_rows(distances_path)

        graph = nx.Graph()
        for row in hap_rows:
            graph.add_node(row["Haplotype"])
        for row in distances:
            graph.add_edge(row["Haplotype1"], row["Haplotype2"], weight=int(row["Distance"]))

        mst = nx.minimum_spanning_tree(graph)
        rows = []
        for source, target, data in mst.edges(data=True):
            rows.append({"Source": source, "Target": target, "Distance": str(data["weight"])})

        with out_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["Source", "Target", "Distance"], delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

        self.logger.info("Built MST with %d edges", len(rows))
        self._mark_global_done("compute_minimum_spanning_tree")
        return out_path

    def draw_network(self, haplotypes_path: Path, samples_path: Path, mst_path: Path) -> Path:
        out_path = self.output_dir / "hap_network.svg"
        if self._should_skip_global("draw_network", [out_path]):
            self.logger.info("Skipping haplotype network drawing")
            return out_path

        try:
            import matplotlib.pyplot as plt
            import networkx as nx
            import numpy as np
            import pandas as pd
        except ImportError as exc:
            raise RuntimeError(
                "Drawing the haplotype network requires matplotlib, networkx, numpy, and pandas"
            ) from exc

        hap_df = pd.read_csv(haplotypes_path, sep="\t")
        sample_df = pd.read_csv(samples_path, sep="\t")
        mst_df = pd.read_csv(mst_path, sep="\t")

        group_col = self.args.group_column
        if group_col not in sample_df.columns:
            raise ValueError(f"Group column '{group_col}' not found in resource metadata")

        regions = list(sample_df[group_col].dropna().astype(str).unique())
        if self.args.group_order:
            requested = [item.strip() for item in self.args.group_order.split(",") if item.strip()]
            regions = requested + [item for item in regions if item not in requested]
        else:
            regions = sorted(regions)

        graph = nx.Graph()
        for _, row in hap_df.iterrows():
            graph.add_node(row["Haplotype"])
        for _, row in mst_df.iterrows():
            graph.add_edge(row["Source"], row["Target"], weight=int(row["Distance"]))

        fig = plt.figure(figsize=(10, 10))
        ax = plt.gca()
        pos = nx.spring_layout(
            graph,
            k=self.args.network_k,
            iterations=self.args.network_iterations,
            seed=self.args.network_seed,
        )

        nx.draw_networkx_edges(graph, pos, width=1.0, edge_color="#D3D3D3", alpha=1.0, ax=ax)

        cmap = plt.get_cmap("Set3")
        region_colors = [cmap(i) for i in np.linspace(0, 1, max(1, len(regions)))]
        color_map = {region: region_colors[idx] for idx, region in enumerate(regions)}

        sample_group_map = dict(zip(sample_df["Sample"], sample_df[group_col].astype(str)))
        hap_sample_map = dict(zip(hap_df["Haplotype"], hap_df["Sample_List"]))
        total_samples = sample_df["Sample"].nunique()

        for _, hap_row in hap_df.iterrows():
            hap = hap_row["Haplotype"]
            sample_names = [item.strip() for item in str(hap_sample_map[hap]).split(",") if item.strip()]
            node_groups = [sample_group_map.get(sample, "Unknown") for sample in sample_names]
            counts = [node_groups.count(region) for region in regions]

            freq = len(sample_names) / max(1, total_samples)
            radius = 0.08 + np.sqrt(freq) * 0.35

            plt.pie(
                counts,
                center=pos[hap],
                radius=radius,
                colors=[color_map[region] for region in regions],
                wedgeprops={"edgecolor": "white", "linewidth": 0.5},
            )

            plt.text(
                pos[hap][0],
                pos[hap][1] + radius + 0.03,
                hap,
                ha="center",
                va="bottom",
                fontsize=12,
                fontweight="bold",
            )

        legend_elements = [
            plt.Line2D([0], [0], marker="o", color="w", label=region, markerfacecolor=color_map[region], markersize=12)
            for region in regions
        ]
        ax.legend(handles=legend_elements, title=group_col, loc="upper right", frameon=True)
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(out_path, format="svg", bbox_inches="tight")
        plt.close(fig)

        self._mark_global_done("draw_network")
        return out_path

    def _read_resource_rows(self, resource_path: Path) -> list[dict[str, str]]:
        if not resource_path.exists():
            raise FileNotFoundError(f"Missing resource file: {resource_path}")
        with resource_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
        if not rows:
            raise ValueError("Resource file is empty")
        if "Sample" not in rows[0]:
            raise ValueError("Resource file must contain a 'Sample' column")
        return rows

    def _read_tsv_rows(self, path: Path) -> list[dict[str, str]]:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return list(reader)

    def _ordered_fieldnames(self, rows: list[dict[str, str]], preferred: list[str]) -> list[str]:
        all_fields = set()
        for row in rows:
            all_fields.update(row.keys())
        ordered = [field for field in preferred if field in all_fields]
        ordered.extend(sorted(field for field in all_fields if field not in ordered))
        return ordered

    def _get_iupac(self, gt_str: str, alt_alleles: list[str]) -> str:
        sep = "|" if "|" in gt_str else "/"
        indices = gt_str.split(sep)
        try:
            a1 = alt_alleles[int(indices[0])] if indices[0] != "." else "."
            a2 = alt_alleles[int(indices[1])] if indices[1] != "." else "."
            return IUPAC.get((a1, a2), "N")
        except Exception:
            return "N"

    def _hamming_distance(self, seq1: str, seq2: str) -> int:
        return sum(base1 != base2 for base1, base2 in zip(seq1, seq2))

    def _ensure_dirs(self) -> dict[str, Path]:
        dirs = {
            "base": self.output_dir,
            "logs": self.output_dir / "logs",
            "state": self.output_dir / "state",
        }
        for path in dirs.values():
            path.mkdir(parents=True, exist_ok=True)
        return dirs

    def _build_logger(self) -> logging.Logger:
        logger = logging.getLogger(f"gwas_step7_{id(self)}")
        logger.setLevel(logging.INFO)
        logger.propagate = False
        if logger.handlers:
            return logger

        formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
        file_handler = logging.FileHandler(self.dirs["logs"] / "pipeline.log")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)
        return logger

    def _load_state(self) -> dict[str, Any]:
        if self.state_path.exists():
            return json.loads(self.state_path.read_text())
        return {"global": {}}

    def _save_state(self) -> None:
        self.state_path.write_text(json.dumps(self.state, indent=2, sort_keys=True))

    def _mark_global_done(self, step: str) -> None:
        self.state["global"][step] = "done"
        self._save_state()

    def _should_skip_global(self, step: str, outputs: list[Path]) -> bool:
        if self.args.force:
            return False
        state_done = self.state["global"].get(step) == "done"
        files_done = all(path.exists() for path in outputs)
        return state_done or files_done


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="GWAS step7 pipeline: haplotype network from VCF and resource metadata."
    )
    parser.add_argument("--vcf", required=True, help="Regional VCF file")
    parser.add_argument("--resource", required=True, help="Resource CSV file")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--group-column", default="Region", help="Resource metadata column used for pie slices")
    parser.add_argument(
        "--group-order",
        help="Optional comma-separated group order, for example India,China,Japan,Russia",
    )
    parser.add_argument("--network-k", type=float, default=3.2, help="Spring layout k parameter")
    parser.add_argument("--network-seed", type=int, default=263, help="Spring layout random seed")
    parser.add_argument("--network-iterations", type=int, default=200, help="Spring layout iterations")
    parser.add_argument("--force", action="store_true", help="Rerun even if outputs already exist")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    pipeline = Step7Pipeline(args)
    pipeline.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
