from __future__ import annotations

import argparse
import csv
import json
import logging
from collections import Counter, defaultdict
from pathlib import Path
import sys
from typing import Any


DEFAULT_COLORS = [
    "#66C2A5",
    "#FC8D62",
    "#8DA0CB",
    "#E78AC3",
    "#A6D854",
    "#FFD92F",
    "#E5C494",
    "#B3B3B3",
    "#1F78B4",
    "#FB9A99",
]


class Step6Pipeline:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.output_dir = Path(args.output_dir).resolve()
        self.dirs = self._ensure_dirs()
        self.logger = self._build_logger()
        self.state_path = self.dirs["state"] / "pipeline_state.json"
        self.state = self._load_state()

    def run(self) -> None:
        self.logger.info("Starting GWAS step6 haplotype analysis")
        hap_csv = self.vcf_to_haplotypes()
        sample_table = self.expand_haplotype_samples(hap_csv)
        group_summary = self.build_group_summary(sample_table)
        self.plot_group_summary(group_summary)
        self.logger.info("Step6 workflow finished successfully")

    def vcf_to_haplotypes(self) -> Path:
        out_path = self.output_dir / "gene.csv"
        if self._should_skip_global("vcf_to_haplotypes", [out_path]):
            self.logger.info("Skipping VCF to haplotype conversion")
            return out_path

        vcf_path = Path(self.args.vcf).resolve()
        if not vcf_path.exists():
            raise FileNotFoundError(f"Missing VCF file: {vcf_path}")

        samples: list[str] = []
        sample_genotypes: dict[str, list[str]] = {}

        with vcf_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    samples = line.strip().split("\t")[9:]
                    sample_genotypes = {sample: [] for sample in samples}
                    continue

                cols = line.strip().split("\t")
                if len(cols) < 10:
                    continue
                ref = cols[3]
                alt_alleles = [ref] + cols[4].split(",")
                genotype_data = cols[9:]

                for idx, gt_field in enumerate(genotype_data):
                    gt = gt_field.split(":", 1)[0]
                    sep = "|" if "|" in gt else "/"
                    indices = gt.split(sep)
                    genotype_str = self._get_allele(indices, alt_alleles)
                    sample_genotypes[samples[idx]].append(genotype_str)

        final_combinations = {}
        for sample in samples:
            final_combinations[sample] = "-".join(sample_genotypes[sample])

        counts = Counter(final_combinations.values())
        sorted_combinations = counts.most_common()

        comb_to_samples: dict[str, list[str]] = defaultdict(list)
        for sample, comb in final_combinations.items():
            comb_to_samples[comb].append(sample)

        total_samples = len(samples)
        rows = []
        for idx, (comb, count) in enumerate(sorted_combinations, start=1):
            rows.append(
                {
                    "Haplotype": f"Hap{idx}",
                    "Genotype_Combination": comb,
                    "Count": str(count),
                    "Frequency": f"{(count / total_samples) * 100:.2f}%",
                    "Sample_List": ",".join(comb_to_samples[comb]),
                }
            )

        with out_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["Haplotype", "Genotype_Combination", "Count", "Frequency", "Sample_List"],
            )
            writer.writeheader()
            writer.writerows(rows)

        self.logger.info("Wrote %d haplotypes ranked by frequency", len(rows))
        self._mark_global_done("vcf_to_haplotypes")
        return out_path

    def expand_haplotype_samples(self, hap_csv: Path) -> Path:
        out_path = self.output_dir / "haplotype_sample_table.tsv"
        if self._should_skip_global("expand_haplotype_samples", [out_path]):
            self.logger.info("Skipping haplotype sample expansion")
            return out_path

        resource_path = Path(self.args.resource).resolve()
        resources = self._read_resource_table(resource_path)

        rows = []
        with hap_csv.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                hap = row["Haplotype"]
                sample_list = [sample.strip() for sample in row["Sample_List"].split(",") if sample.strip()]
                for sample in sample_list:
                    meta = resources.get(sample, {})
                    out_row = {"Sample": sample, "Haplotype": hap}
                    for key, value in meta.items():
                        out_row[key] = value
                    rows.append(out_row)

        fieldnames = self._ordered_fieldnames(rows, ["Sample", "Haplotype"])
        with out_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

        self.logger.info("Expanded haplotypes to %d sample-level rows", len(rows))
        self._mark_global_done("expand_haplotype_samples")
        return out_path

    def build_group_summary(self, sample_table: Path) -> Path:
        out_path = self.output_dir / "haplotype_group_summary.tsv"
        if self._should_skip_global("build_group_summary", [out_path]):
            self.logger.info("Skipping group summary")
            return out_path

        with sample_table.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)

        group_col = self.args.group_column
        if not rows:
            raise ValueError("No sample rows available for group summary")
        if group_col not in rows[0]:
            raise ValueError(f"Group column '{group_col}' not found in resource metadata")

        group_hap_counts: dict[str, Counter[str]] = defaultdict(Counter)
        hap_order = []
        seen_haps = set()

        for row in rows:
            group = row[group_col]
            hap = row["Haplotype"]
            group_hap_counts[group][hap] += 1
            if hap not in seen_haps:
                hap_order.append(hap)
                seen_haps.add(hap)

        summary_rows = []
        for group in self._group_order(group_hap_counts.keys()):
            total = sum(group_hap_counts[group].values())
            for hap in hap_order:
                count = group_hap_counts[group].get(hap, 0)
                proportion = (count / total * 100.0) if total else 0.0
                summary_rows.append(
                    {
                        group_col: group,
                        "Haplotype": hap,
                        "Count": str(count),
                        "Proportion": f"{proportion:.6f}",
                    }
                )

        with out_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[group_col, "Haplotype", "Count", "Proportion"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(summary_rows)

        self.logger.info("Wrote haplotype summary across %d groups", len(group_hap_counts))
        self._mark_global_done("build_group_summary")
        return out_path

    def plot_group_summary(self, summary_path: Path) -> Path:
        out_path = self.dirs["plots"] / "haplotype_by_group.svg"
        if self._should_skip_global("plot_group_summary", [out_path]):
            self.logger.info("Skipping haplotype group plot")
            return out_path

        with summary_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)

        if not rows:
            raise ValueError("No summary rows available for plotting")

        group_col = self.args.group_column
        groups = []
        haplotypes = []
        seen_groups = set()
        seen_haps = set()
        for row in rows:
            if row[group_col] not in seen_groups:
                groups.append(row[group_col])
                seen_groups.add(row[group_col])
            if row["Haplotype"] not in seen_haps:
                haplotypes.append(row["Haplotype"])
                seen_haps.add(row["Haplotype"])

        summary_map: dict[tuple[str, str], float] = {}
        for row in rows:
            summary_map[(row[group_col], row["Haplotype"])] = float(row["Proportion"])

        width = 1100
        height = 700
        left = 100
        right = 160
        top = 80
        bottom = 90
        plot_width = width - left - right
        plot_height = height - top - bottom

        elements = ['<rect width="100%" height="100%" fill="white" />']
        elements.append(f'<text x="{left}" y="34" font-size="28" font-weight="700">Haplotype Composition by {self._escape(group_col)}</text>')
        elements.append(f'<text x="{left}" y="60" font-size="15" fill="#666">Proportional stacked bar chart derived from sample metadata.</text>')

        for pct in [0, 25, 50, 75, 100]:
            y = top + plot_height - (pct / 100.0) * plot_height
            elements.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" stroke="#ebebeb" />')
            elements.append(f'<text x="{left - 10}" y="{y + 5:.2f}" text-anchor="end" font-size="12">{pct}%</text>')

        elements.append(f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#333" />')
        elements.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#333" />')
        elements.append(f'<text x="{left + plot_width / 2:.2f}" y="{height - 18}" text-anchor="middle" font-size="16">{self._escape(group_col)}</text>')
        elements.append(f'<text x="26" y="{top + plot_height / 2:.2f}" text-anchor="middle" font-size="16" transform="rotate(-90 26 {top + plot_height / 2:.2f})">Proportion (%)</text>')

        slot = plot_width / max(1, len(groups))
        bar_width = slot * 0.56
        color_map = {hap: DEFAULT_COLORS[idx % len(DEFAULT_COLORS)] for idx, hap in enumerate(haplotypes)}

        for idx, group in enumerate(groups):
            x = left + idx * slot + (slot - bar_width) / 2
            cumulative = 0.0
            for hap in haplotypes:
                proportion = summary_map.get((group, hap), 0.0)
                seg_height = (proportion / 100.0) * plot_height
                y = top + plot_height - cumulative - seg_height
                elements.append(
                    f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_width:.2f}" height="{seg_height:.2f}" fill="{color_map[hap]}" />'
                )
                if proportion >= 7:
                    elements.append(
                        f'<text x="{x + bar_width / 2:.2f}" y="{y + seg_height / 2 + 4:.2f}" '
                        f'text-anchor="middle" font-size="12" fill="white">{proportion:.1f}%</text>'
                    )
                cumulative += seg_height
            elements.append(f'<text x="{x + bar_width / 2:.2f}" y="{top + plot_height + 24:.2f}" text-anchor="middle" font-size="13">{self._escape(group)}</text>')

        legend_x = left + plot_width + 30
        legend_y = top + 10
        for idx, hap in enumerate(haplotypes):
            y = legend_y + idx * 24
            elements.append(f'<rect x="{legend_x}" y="{y}" width="16" height="16" fill="{color_map[hap]}" rx="3" />')
            elements.append(f'<text x="{legend_x + 24}" y="{y + 13}" font-size="13">{self._escape(hap)}</text>')

        svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">
{''.join(elements)}
</svg>
"""
        out_path.write_text(svg, encoding="utf-8")
        self._mark_global_done("plot_group_summary")
        return out_path

    def _read_resource_table(self, resource_path: Path) -> dict[str, dict[str, str]]:
        if not resource_path.exists():
            raise FileNotFoundError(f"Missing resource file: {resource_path}")

        with resource_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
        if not rows:
            raise ValueError("Resource file is empty")
        if "Sample" not in rows[0]:
            raise ValueError("Resource file must contain a 'Sample' column")
        return {row["Sample"]: row for row in rows}

    def _group_order(self, groups: Any) -> list[str]:
        discovered = list(groups)
        if self.args.group_order:
            requested = [item.strip() for item in self.args.group_order.split(",") if item.strip()]
            remaining = [group for group in discovered if group not in requested]
            return requested + remaining
        return sorted(discovered)

    def _ordered_fieldnames(self, rows: list[dict[str, str]], preferred: list[str]) -> list[str]:
        all_fields = set()
        for row in rows:
            all_fields.update(row.keys())
        ordered = [field for field in preferred if field in all_fields]
        ordered.extend(sorted(field for field in all_fields if field not in ordered))
        return ordered

    def _get_allele(self, gt_indices: list[str], alt_alleles: list[str]) -> str:
        bases = []
        for idx in gt_indices:
            if idx == ".":
                bases.append(".")
            else:
                try:
                    bases.append(alt_alleles[int(idx)])
                except (ValueError, IndexError):
                    bases.append("N")
        return "/".join(bases)

    def _ensure_dirs(self) -> dict[str, Path]:
        dirs = {
            "base": self.output_dir,
            "plots": self.output_dir / "plots",
            "logs": self.output_dir / "logs",
            "state": self.output_dir / "state",
        }
        for path in dirs.values():
            path.mkdir(parents=True, exist_ok=True)
        return dirs

    def _build_logger(self) -> logging.Logger:
        logger = logging.getLogger(f"gwas_step6_{id(self)}")
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

    @staticmethod
    def _escape(text: str) -> str:
        return (
            str(text)
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace('"', "&quot;")
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="GWAS step6 pipeline: haplotype annotation and metadata-based composition plot."
    )
    parser.add_argument("--vcf", required=True, help="Regional VCF file")
    parser.add_argument("--resource", required=True, help="Sample resource metadata CSV")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--group-column", default="Region", help="Metadata column used for grouping")
    parser.add_argument(
        "--group-order",
        help="Optional comma-separated order for group display, for example India,China,Japan,Russia",
    )
    parser.add_argument("--force", action="store_true", help="Rerun even if outputs already exist")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    pipeline = Step6Pipeline(args)
    pipeline.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
