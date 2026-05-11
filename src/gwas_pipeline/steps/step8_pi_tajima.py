from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
from pathlib import Path
import re
import shlex
import subprocess
import sys
from typing import Any, Optional


REGION_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")


class Step8Pipeline:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.output_dir = Path(args.output_dir).resolve()
        self.dirs = self._ensure_dirs()
        self.logger = self._build_logger()
        self.state_path = self.dirs["state"] / "pipeline_state.json"
        self.state = self._load_state()
        self.region = self._parse_region(args.region)

    def run(self) -> None:
        self.logger.info("Starting GWAS step8 pi and Tajima's D workflow")
        region_vcf = self.prepare_region_vcf()
        groups = self.read_groups()
        group_vcfs = self.split_group_vcfs(region_vcf, groups)
        self.compute_population_stats(group_vcfs)
        pi_merged = self.merge_pi_results(group_vcfs)
        tajima_merged = self.merge_tajima_results(group_vcfs)
        self.plot_results(pi_merged, tajima_merged)
        self.logger.info("Step8 workflow finished successfully")

    def prepare_region_vcf(self) -> Path:
        source_vcf = Path(self.args.vcf).resolve()
        if not source_vcf.exists():
            raise FileNotFoundError(f"Missing VCF file: {source_vcf}")

        out_path = self.output_dir / "gene_region.vcf.gz"
        if self._should_skip_global("prepare_region_vcf", [out_path]):
            self.logger.info("Skipping regional VCF extraction")
            return out_path

        cmd = (
            f"bcftools view {self._q(source_vcf)} "
            f"-r {shlex.quote(self.args.region)} "
            f"-Oz -o {self._q(out_path)}"
        )
        self._run_cmd(cmd)
        self._run_cmd(f"bcftools index -t {self._q(out_path)}")
        out_path = self._ensure_vcftools_compatible_vcf(out_path, "prepare_region_vcf_v42")
        self._mark_global_done("prepare_region_vcf")
        return out_path

    def read_groups(self) -> dict[str, list[str]]:
        class_path = Path(self.args.class_file).resolve()
        if not class_path.exists():
            raise FileNotFoundError(f"Missing class file: {class_path}")

        with class_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
        if not rows:
            raise ValueError("Class file is empty")

        sample_col = self.args.sample_column
        group_col = self.args.group_column
        if sample_col not in rows[0]:
            raise ValueError(f"Sample column '{sample_col}' not found")
        if group_col not in rows[0]:
            raise ValueError(f"Group column '{group_col}' not found")

        groups: dict[str, list[str]] = {}
        for row in rows:
            groups.setdefault(row[group_col], []).append(row[sample_col])
        return groups

    def split_group_vcfs(self, region_vcf: Path, groups: dict[str, list[str]]) -> dict[str, Path]:
        outputs = {}
        for group_name, samples in sorted(groups.items()):
            sample_file = self.dirs["groups"] / f"{group_name}_samples.txt"
            group_vcf = self.dirs["groups"] / f"{group_name}.vcf.gz"
            outputs[group_name] = group_vcf

            if not self._should_skip_global(f"write_samples_{group_name}", [sample_file]):
                sample_file.write_text("\n".join(samples) + "\n", encoding="utf-8")
                self._mark_global_done(f"write_samples_{group_name}")

            if self._should_skip_global(f"split_vcf_{group_name}", [group_vcf]):
                self.logger.info("Skipping group VCF split for %s", group_name)
                continue

            cmd = (
                f"bcftools view -S {self._q(sample_file)} "
                f"{self._q(region_vcf)} "
                f"-Oz -o {self._q(group_vcf)}"
            )
            self._run_cmd(cmd)
            self._run_cmd(f"bcftools index -t {self._q(group_vcf)}")
            compatible_vcf = self._ensure_vcftools_compatible_vcf(
                group_vcf, f"split_vcf_{group_name}_v42"
            )
            outputs[group_name] = compatible_vcf
            self._mark_global_done(f"split_vcf_{group_name}")
        return outputs

    def compute_population_stats(self, group_vcfs: dict[str, Path]) -> None:
        chrom = self.region["chrom"]
        start = self.region["start"]
        end = self.region["end"]

        for group_name, group_vcf in sorted(group_vcfs.items()):
            pi_prefix = self.dirs["stats"] / f"{group_name}_pi"
            tajima_prefix = self.dirs["stats"] / f"{group_name}_tajima"
            pi_file = Path(f"{pi_prefix}.windowed.pi")
            tajima_file = Path(f"{tajima_prefix}.Tajima.D")

            if not self._should_skip_global(f"pi_{group_name}", [pi_file]):
                cmd = (
                    f"vcftools --gzvcf {self._q(group_vcf)} "
                    f"--chr {shlex.quote(chrom)} "
                    f"--from-bp {start} "
                    f"--to-bp {end} "
                    f"--window-pi {self.args.window_size} "
                    f"--window-pi-step {self.args.step_size} "
                    f"--out {self._q(pi_prefix)}"
                )
                self._run_cmd(cmd)
                self._mark_global_done(f"pi_{group_name}")

            if not self._should_skip_global(f"tajima_{group_name}", [tajima_file]):
                cmd = (
                    f"vcftools --gzvcf {self._q(group_vcf)} "
                    f"--chr {shlex.quote(chrom)} "
                    f"--from-bp {start} "
                    f"--to-bp {end} "
                    f"--TajimaD {self.args.window_size} "
                    f"--out {self._q(tajima_prefix)}"
                )
                self._run_cmd(cmd)
                self._mark_global_done(f"tajima_{group_name}")

    def merge_pi_results(self, group_vcfs: dict[str, Path]) -> Path:
        out_path = self.dirs["merged"] / "pi_merged.tsv"
        if self._should_skip_global("merge_pi_results", [out_path]):
            self.logger.info("Skipping merged pi table")
            return out_path

        rows = []
        for group_name in sorted(group_vcfs):
            pi_file = self.dirs["stats"] / f"{group_name}_pi.windowed.pi"
            if not pi_file.exists():
                continue
            with pi_file.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    pi_value = row.get("PI", "")
                    if pi_value in {"", "NA", "nan"}:
                        continue
                    try:
                        pi_numeric = float(pi_value)
                    except ValueError:
                        continue
                    if pi_numeric < 0:
                        continue
                    rows.append(
                        {
                            "Group": group_name,
                            "CHROM": row["CHROM"],
                            "BIN_START": row["BIN_START"],
                            "BIN_END": row.get("BIN_END", str(int(row["BIN_START"]) + self.args.window_size)),
                            "Position": str(int(row["BIN_START"]) - 1),
                            "PI": f"{pi_numeric:.10g}",
                        }
                    )

        self._write_tsv(out_path, ["Group", "CHROM", "BIN_START", "BIN_END", "Position", "PI"], rows)
        self._mark_global_done("merge_pi_results")
        return out_path

    def merge_tajima_results(self, group_vcfs: dict[str, Path]) -> Path:
        out_path = self.dirs["merged"] / "tajima_merged.tsv"
        if self._should_skip_global("merge_tajima_results", [out_path]):
            self.logger.info("Skipping merged Tajima table")
            return out_path

        rows = []
        for group_name in sorted(group_vcfs):
            tajima_file = self.dirs["stats"] / f"{group_name}_tajima.Tajima.D"
            if not tajima_file.exists():
                continue
            with tajima_file.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    tajima_value = row.get("TajimaD", "")
                    if tajima_value in {"", "NA", "nan"}:
                        continue
                    try:
                        tajima_numeric = float(tajima_value)
                    except ValueError:
                        continue
                    rows.append(
                        {
                            "Group": group_name,
                            "CHROM": row["CHROM"],
                            "BIN_START": row["BIN_START"],
                            "BIN_END": row.get("BIN_END", str(int(row["BIN_START"]) + self.args.window_size)),
                            "Position": row["BIN_START"],
                            "TajimaD": f"{tajima_numeric:.10g}",
                        }
                    )

        self._write_tsv(out_path, ["Group", "CHROM", "BIN_START", "BIN_END", "Position", "TajimaD"], rows)
        self._mark_global_done("merge_tajima_results")
        return out_path

    def plot_results(self, pi_merged: Path, tajima_merged: Path) -> Path:
        out_path = self.dirs["plots"] / "pi_tajima.svg"
        if self._should_skip_global("plot_results", [out_path]):
            self.logger.info("Skipping pi/Tajima plot")
            return out_path

        pi_rows = self._read_tsv_rows(pi_merged)
        tajima_rows = self._read_tsv_rows(tajima_merged)

        pi_rows = [row for row in pi_rows if row["PI"] not in {"", "NA"}]
        tajima_rows = [row for row in tajima_rows if row["TajimaD"] not in {"", "NA"}]

        sig_regions = self._detect_sig_regions(tajima_rows)
        svg = self._build_plot_svg(pi_rows, tajima_rows, sig_regions)
        out_path.write_text(svg, encoding="utf-8")
        self._mark_global_done("plot_results")
        return out_path

    def _detect_sig_regions(self, tajima_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
        rows = []
        grouped: dict[str, list[dict[str, Any]]] = {}
        for row in tajima_rows:
            grouped.setdefault(row["Group"], []).append(
                {
                    "Position": int(row["Position"]),
                    "TajimaD": float(row["TajimaD"]),
                }
            )

        for group_name, items in grouped.items():
            hits = [item for item in sorted(items, key=lambda x: x["Position"]) if item["TajimaD"] < -2]
            if not hits:
                continue
            block_start = hits[0]["Position"]
            prev = hits[0]["Position"]
            for item in hits[1:]:
                if item["Position"] - prev > (2 * self.args.window_size):
                    rows.append({"Group": group_name, "start": block_start, "end": prev + self.args.window_size})
                    block_start = item["Position"]
                prev = item["Position"]
            rows.append({"Group": group_name, "start": block_start, "end": prev + self.args.window_size})
        return rows

    def _build_plot_svg(
        self,
        pi_rows: list[dict[str, str]],
        tajima_rows: list[dict[str, str]],
        sig_regions: list[dict[str, Any]],
    ) -> str:
        width = 1200
        height = 820
        left = 90
        right = 40
        top = 70
        gap = 50
        bottom = 70
        panel_height = (height - top - bottom - gap) / 2
        plot_width = width - left - right
        xmin = self.region["start"]
        xmax = self.region["end"]

        colors = {
            "Cultivar": "#E67E22",
            "Landrace": "#7A9AAB",
            "Wild": "#2E5984",
        }

        ymax_pi = max((float(row["PI"]) for row in pi_rows), default=0.01) * 1.2
        ymin_d = min((float(row["TajimaD"]) for row in tajima_rows), default=-3.0) - 0.5
        ymax_d = max(2.0, max((float(row["TajimaD"]) for row in tajima_rows), default=2.0) + 0.5)

        elements = ['<rect width="100%" height="100%" fill="white" />']
        elements.append(f'<text x="{left}" y="32" font-size="26" font-weight="700">Pi and Tajima&apos;s D in target region</text>')

        elements.extend(self._draw_panel_frame(left, top, plot_width, panel_height, "Nucleotide diversity", "π"))
        elements.extend(self._draw_panel_frame(left, top + panel_height + gap, plot_width, panel_height, "", "Tajima&apos;s D"))

        elements.extend(self._draw_x_ticks(left, top + panel_height + gap + panel_height, plot_width, xmin, xmax))
        elements.append(
            f'<text x="{left + plot_width / 2:.2f}" y="{height - 18}" text-anchor="middle" font-size="16">Genomic position</text>'
        )

        for group_name in self._ordered_groups(pi_rows, tajima_rows):
            pi_points = [(int(row["Position"]), float(row["PI"])) for row in pi_rows if row["Group"] == group_name]
            tajima_points = [(int(row["Position"]), float(row["TajimaD"])) for row in tajima_rows if row["Group"] == group_name]
            color = colors.get(group_name, "#555555")
            elements.append(self._polyline(pi_points, left, top, plot_width, panel_height, xmin, xmax, 0, ymax_pi, color))
            elements.append(self._polyline(tajima_points, left, top + panel_height + gap, plot_width, panel_height, xmin, xmax, ymin_d, ymax_d, color))

        for region in sig_regions:
            color = colors.get(region["Group"], "#cccccc")
            x1 = left + ((region["start"] - xmin) / max(1, xmax - xmin)) * plot_width
            x2 = left + ((region["end"] - xmin) / max(1, xmax - xmin)) * plot_width
            y = top + panel_height + gap
            elements.append(
                f'<rect x="{x1:.2f}" y="{y:.2f}" width="{max(1.0, x2 - x1):.2f}" height="{panel_height:.2f}" fill="{color}" fill-opacity="0.18" />'
            )

        threshold_y = top + panel_height + gap + panel_height - ((-2 - ymin_d) / max(1e-9, (ymax_d - ymin_d))) * panel_height
        elements.append(
            f'<line x1="{left}" y1="{threshold_y:.2f}" x2="{left + plot_width}" y2="{threshold_y:.2f}" stroke="red" stroke-dasharray="6,4" />'
        )

        legend_x = width - right - 130
        legend_y = top - 5
        for idx, group_name in enumerate(self._ordered_groups(pi_rows, tajima_rows)):
            y = legend_y + idx * 24
            color = colors.get(group_name, "#555555")
            elements.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x + 18}" y2="{y}" stroke="{color}" stroke-width="3" />')
            elements.append(f'<text x="{legend_x + 26}" y="{y + 4}" font-size="13">{self._escape(group_name)}</text>')

        return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">
{''.join(elements)}
</svg>
"""

    def _draw_panel_frame(self, left: int, top: float, width: int, height: float, title: str, ylabel: str) -> list[str]:
        parts = [
            f'<line x1="{left}" y1="{top + height}" x2="{left + width}" y2="{top + height}" stroke="#333" />',
            f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + height}" stroke="#333" />',
            f'<text x="26" y="{top + height / 2:.2f}" text-anchor="middle" font-size="16" transform="rotate(-90 26 {top + height / 2:.2f})">{ylabel}</text>',
        ]
        if title:
            parts.append(f'<text x="{left + width / 2:.2f}" y="{top - 12}" text-anchor="middle" font-size="18" font-weight="700">{title}</text>')
        return parts

    def _draw_x_ticks(self, left: int, base_y: float, width: int, xmin: int, xmax: int) -> list[str]:
        parts = []
        for frac in [0.0, 0.25, 0.5, 0.75, 1.0]:
            pos = int(xmin + frac * (xmax - xmin))
            x = left + frac * width
            parts.append(f'<line x1="{x:.2f}" y1="{base_y:.2f}" x2="{x:.2f}" y2="{base_y + 6:.2f}" stroke="#333" />')
            parts.append(f'<text x="{x:.2f}" y="{base_y + 24:.2f}" text-anchor="middle" font-size="12">{pos}</text>')
        return parts

    def _polyline(
        self,
        points: list[tuple[int, float]],
        left: int,
        top: float,
        width: int,
        height: float,
        xmin: int,
        xmax: int,
        ymin: float,
        ymax: float,
        color: str,
    ) -> str:
        if not points:
            return ""
        coords = []
        denom_x = max(1, xmax - xmin)
        denom_y = max(1e-9, ymax - ymin)
        for xval, yval in sorted(points):
            x = left + ((xval - xmin) / denom_x) * width
            y = top + height - ((yval - ymin) / denom_y) * height
            coords.append(f"{x:.2f},{y:.2f}")
        return f'<polyline fill="none" stroke="{color}" stroke-width="2" points="{" ".join(coords)}" />'

    def _ordered_groups(self, pi_rows: list[dict[str, str]], tajima_rows: list[dict[str, str]]) -> list[str]:
        groups = []
        seen = set()
        for row in pi_rows + tajima_rows:
            group = row["Group"]
            if group not in seen:
                groups.append(group)
                seen.add(group)
        preferred = ["Cultivar", "Landrace", "Wild"]
        ordered = [g for g in preferred if g in seen]
        ordered.extend(g for g in groups if g not in ordered)
        return ordered

    def _parse_region(self, region: str) -> dict[str, Any]:
        match = REGION_RE.match(region)
        if not match:
            raise ValueError(f"Invalid region format: {region}")
        return {
            "chrom": match.group("chrom"),
            "start": int(match.group("start")),
            "end": int(match.group("end")),
        }

    def _ensure_dirs(self) -> dict[str, Path]:
        dirs = {
            "base": self.output_dir,
            "groups": self.output_dir / "groups",
            "stats": self.output_dir / "stats",
            "merged": self.output_dir / "merged",
            "plots": self.output_dir / "plots",
            "logs": self.output_dir / "logs",
            "state": self.output_dir / "state",
        }
        for path in dirs.values():
            path.mkdir(parents=True, exist_ok=True)
        return dirs

    def _build_logger(self) -> logging.Logger:
        logger = logging.getLogger(f"gwas_step8_{id(self)}")
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

    def _ensure_vcftools_compatible_vcf(self, source_vcf: Path, state_key: str) -> Path:
        compatible_vcf = source_vcf.with_name(source_vcf.name.replace(".vcf.gz", ".v42.vcf.gz"))
        compatible_tbi = Path(f"{compatible_vcf}.tbi")

        if self._should_skip_global(state_key, [compatible_vcf, compatible_tbi]):
            self.logger.info("Skipping VCFv4.2 compatibility conversion for %s", source_vcf.name)
            return compatible_vcf

        cmd = (
            f"bcftools view {self._q(source_vcf)} | "
            f"sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' | "
            f"bgzip -c > {self._q(compatible_vcf)}"
        )
        self._run_cmd(cmd)
        self._run_cmd(f"tabix -p vcf {self._q(compatible_vcf)}")
        self._mark_global_done(state_key)
        return compatible_vcf

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

    def _run_cmd(self, cmd: str) -> None:
        self.logger.info("Running command: %s", cmd)
        result = subprocess.run(
            cmd,
            shell=True,
            executable="/bin/bash",
            text=True,
            capture_output=True,
        )
        if result.stdout:
            self.logger.info(result.stdout.rstrip())
        if result.stderr:
            self.logger.info(result.stderr.rstrip())
        if result.returncode != 0:
            self.logger.error("Command failed with exit code %s", result.returncode)
            raise RuntimeError(f"Command failed: {cmd}")

    def _write_tsv(self, path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def _read_tsv_rows(self, path: Path) -> list[dict[str, str]]:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return list(reader)

    @staticmethod
    def _q(path: Path) -> str:
        return shlex.quote(str(path))

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
        description="GWAS step8 pipeline: pi and Tajima's D analysis by population group."
    )
    parser.add_argument("--vcf", required=True, help="Input VCF.gz file")
    parser.add_argument("--region", required=True, help="Target region such as Chr3:66400000-66411250")
    parser.add_argument("--class-file", required=True, help="Classification CSV file")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--sample-column", default="Sample", help="Sample column in class file")
    parser.add_argument("--group-column", default="Group", help="Group column in class file")
    parser.add_argument("--window-size", type=int, default=100, help="Window size for pi and Tajima's D")
    parser.add_argument("--step-size", type=int, default=100, help="Step size for pi windows")
    parser.add_argument("--force", action="store_true", help="Rerun even if outputs already exist")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    pipeline = Step8Pipeline(args)
    pipeline.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
