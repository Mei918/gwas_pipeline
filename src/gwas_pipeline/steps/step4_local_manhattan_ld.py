from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path
import shlex
import subprocess
import sys
from typing import Any


class Step4Pipeline:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.output_dir = Path(args.output_dir).resolve()
        self.dirs = self._ensure_dirs()
        self.logger = self._build_logger()
        self.state_path = self.dirs["state"] / "pipeline_state.json"
        self.state = self._load_state()

    def run(self) -> None:
        self.logger.info("Starting GWAS step4 regional visualization workflow")
        region_file = self.extract_region()
        rows = self.load_region_rows(region_file)
        self.write_basic_local_plot(rows)

        if self.args.gene_bed:
            genes = self.load_genes()
            self.write_gene_local_plot(rows, genes)

        if self.args.vcf and self.args.ldblockshow:
            self.run_ldblockshow()

        self.logger.info("Step4 workflow finished successfully")

    def extract_region(self) -> Path:
        out_path = self.output_dir / "region_gwas.tsv"
        if self._should_skip_global("extract_region", [out_path]):
            self.logger.info("Skipping regional GWAS extraction")
            return out_path

        gwas_path = Path(self.args.gwas).resolve()
        with gwas_path.open("r", encoding="utf-8") as handle:
            header_line = handle.readline().strip()
            if not header_line:
                raise ValueError("GWAS file is empty")
            headers = header_line.split()

            for column in [self.args.chr_column, self.args.pos_column, self.args.p_column]:
                if column not in headers:
                    raise ValueError(f"Column '{column}' not found in GWAS file")

            chr_idx = headers.index(self.args.chr_column)
            pos_idx = headers.index(self.args.pos_column)

            out_lines = ["\t".join(headers)]
            kept = 0
            for raw in handle:
                fields = raw.strip().split()
                if len(fields) != len(headers):
                    continue
                chrom = self._normalize_chromosome(fields[chr_idx])
                pos = int(fields[pos_idx])
                if chrom == self._target_chromosome() and self.args.start <= pos <= self.args.end:
                    out_lines.append("\t".join(fields))
                    kept += 1

        out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
        self.logger.info("Extracted %d rows in target region", kept)
        self._mark_global_done("extract_region")
        return out_path

    def load_region_rows(self, region_file: Path) -> list[dict[str, Any]]:
        rows = []
        with region_file.open("r", encoding="utf-8") as handle:
            header = handle.readline().strip().split("\t")
            if len(header) == 1:
                header = handle.readline().strip().split()
            if not header:
                raise ValueError("Region GWAS file has no header")

        with region_file.open("r", encoding="utf-8") as handle:
            header_line = handle.readline().strip()
            headers = header_line.split("\t")
            if len(headers) == 1:
                headers = header_line.split()
            chr_idx = headers.index(self.args.chr_column)
            pos_idx = headers.index(self.args.pos_column)
            p_idx = headers.index(self.args.p_column)
            for raw in handle:
                fields = raw.strip().split("\t")
                if len(fields) != len(headers):
                    fields = raw.strip().split()
                if len(fields) != len(headers):
                    continue
                p_value = max(float(fields[p_idx]), 1e-300)
                rows.append(
                    {
                        "chr": self._normalize_chromosome(fields[chr_idx]),
                        "bp": int(fields[pos_idx]),
                        "p": p_value,
                        "logp": -math.log10(p_value),
                    }
                )
        if not rows:
            self.logger.warning("No GWAS rows found in target region")
        return rows

    def load_genes(self) -> list[dict[str, Any]]:
        gene_bed = Path(self.args.gene_bed).resolve()
        genes = []
        with gene_bed.open("r", encoding="utf-8") as handle:
            for raw in handle:
                fields = raw.strip().split("\t")
                if len(fields) < 6:
                    continue
                chrom = self._normalize_chromosome(fields[0])
                start1 = int(fields[1]) + 1
                end1 = int(fields[2])
                if chrom != self._target_chromosome():
                    continue
                if end1 < self.args.start or start1 > self.args.end:
                    continue
                genes.append(
                    {
                        "chr": chrom,
                        "start": start1,
                        "end": end1,
                        "gene_id": fields[3],
                        "strand": fields[5],
                    }
                )
        return self._assign_gene_rows(sorted(genes, key=lambda item: (item["start"], item["end"])))

    def write_basic_local_plot(self, rows: list[dict[str, Any]]) -> Path:
        out_path = self.dirs["plots"] / "local_manhattan.svg"
        if self._should_skip_global("plot_basic", [out_path]):
            self.logger.info("Skipping basic local Manhattan plot")
            return out_path
        svg = self._build_local_plot_svg(rows, genes=None)
        out_path.write_text(svg, encoding="utf-8")
        self._mark_global_done("plot_basic")
        return out_path

    def write_gene_local_plot(self, rows: list[dict[str, Any]], genes: list[dict[str, Any]]) -> Path:
        out_path = self.dirs["plots"] / "local_manhattan_with_genes.svg"
        if self._should_skip_global("plot_with_genes", [out_path]):
            self.logger.info("Skipping local Manhattan plot with genes")
            return out_path
        svg = self._build_local_plot_svg(rows, genes=genes)
        out_path.write_text(svg, encoding="utf-8")
        self._mark_global_done("plot_with_genes")
        return out_path

    def run_ldblockshow(self) -> Path:
        output_prefix = self.dirs["ld"] / "test_output"
        marker = Path(f"{output_prefix}.site.gz")
        if self._should_skip_global("run_ldblockshow", [marker]):
            self.logger.info("Skipping LDBlockShow")
            return output_prefix

        vcf = Path(self.args.vcf).resolve()
        exe = Path(self.args.ldblockshow).resolve()
        region = f"{self._target_chromosome()}:{self.args.start}-{self.args.end}"
        cmd = (
            f"{self._q(exe)} "
            f"-InVCF {self._q(vcf)} "
            f"-OutPut {self._q(output_prefix)} "
            f"-Region {shlex.quote(region)}"
        )
        self._run_cmd(cmd)
        self._mark_global_done("run_ldblockshow")
        return output_prefix

    def _build_local_plot_svg(self, rows: list[dict[str, Any]], genes: list[dict[str, Any]] | None) -> str:
        width = 1200
        height = 720 if genes else 560
        left = 90
        right = 35
        top = 75
        bottom = 90 if genes else 70
        plot_width = width - left - right
        plot_height = height - top - bottom

        max_logp = max((row["logp"] for row in rows), default=1.0)
        threshold = -math.log10(self.args.threshold_line)
        max_y = max(max_logp, threshold) * (1.08 if rows else 1.0)
        if max_y < 1.0:
            max_y = 1.0

        background = ['<rect width="100%" height="100%" fill="white" />']
        grid = []
        for frac in [0.25, 0.5, 0.75, 1.0]:
            y = top + plot_height - frac * plot_height
            value = frac * max_y
            grid.append(f'<line x1="{left}" y1="{y:.2f}" x2="{width - right}" y2="{y:.2f}" stroke="#ececec" />')
            grid.append(f'<text x="{left - 10}" y="{y + 5:.2f}" text-anchor="end" font-size="12">{value:.1f}</text>')

        axis = [
            f'<line x1="{left}" y1="{top + plot_height}" x2="{width - right}" y2="{top + plot_height}" stroke="#333" />',
            f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#333" />',
            f'<text x="{width / 2:.0f}" y="34" text-anchor="middle" font-size="26">Local Manhattan Plot</text>',
            f'<text x="{width / 2:.0f}" y="{height - 18}" text-anchor="middle" font-size="16">Position (bp)</text>',
            f'<text x="24" y="{height / 2:.0f}" text-anchor="middle" font-size="16" transform="rotate(-90 24 {height / 2:.0f})">-log10(P)</text>',
        ]

        x_ticks = []
        for frac in [0.0, 0.25, 0.5, 0.75, 1.0]:
            pos = int(self.args.start + frac * (self.args.end - self.args.start))
            x = left + frac * plot_width
            x_ticks.append(f'<line x1="{x:.2f}" y1="{top + plot_height}" x2="{x:.2f}" y2="{top + plot_height + 6}" stroke="#333" />')
            x_ticks.append(f'<text x="{x:.2f}" y="{top + plot_height + 24}" text-anchor="middle" font-size="12">{pos}</text>')

        threshold_y = top + plot_height - (threshold / max_y) * plot_height
        threshold_line = f'<line x1="{left}" y1="{threshold_y:.2f}" x2="{width - right}" y2="{threshold_y:.2f}" stroke="red" stroke-dasharray="6,4" />'

        points = []
        for row in rows:
            frac = (row["bp"] - self.args.start) / max(1, self.args.end - self.args.start)
            x = left + frac * plot_width
            y = top + plot_height - (row["logp"] / max_y) * plot_height
            points.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.0" fill="#3E0A52" fill-opacity="0.82" />')

        gene_shapes = []
        if genes:
            gene_base = top + plot_height + 35
            row_step = 28
            for gene in genes:
                y = gene_base + (gene["row"] - 1) * row_step
                start_x = left + ((gene["start"] - self.args.start) / max(1, self.args.end - self.args.start)) * plot_width
                end_x = left + ((gene["end"] - self.args.start) / max(1, self.args.end - self.args.start)) * plot_width
                start_x = max(left, min(width - right, start_x))
                end_x = max(left, min(width - right, end_x))
                if gene["strand"] == "+":
                    x1, x2 = start_x, end_x
                else:
                    x1, x2 = end_x, start_x
                gene_shapes.append(
                    f'<line x1="{x1:.2f}" y1="{y:.2f}" x2="{x2:.2f}" y2="{y:.2f}" stroke="black" stroke-width="2" />'
                )
                gene_shapes.append(self._arrow_head(x2, y, gene["strand"]))
                mid = (start_x + end_x) / 2
                gene_shapes.append(
                    f'<text x="{mid:.2f}" y="{y + 16:.2f}" text-anchor="middle" font-size="12" fill="red">{self._escape(gene["gene_id"])}</text>'
                )

        return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">
{''.join(background)}
{''.join(grid)}
{''.join(axis)}
{''.join(x_ticks)}
{threshold_line}
{''.join(points)}
{''.join(gene_shapes)}
</svg>
"""

    def _assign_gene_rows(self, genes: list[dict[str, Any]]) -> list[dict[str, Any]]:
        row_ends: list[int] = []
        for gene in genes:
            placed = False
            for idx, row_end in enumerate(row_ends):
                if gene["start"] > row_end:
                    gene["row"] = idx + 1
                    row_ends[idx] = gene["end"]
                    placed = True
                    break
            if not placed:
                row_ends.append(gene["end"])
                gene["row"] = len(row_ends)
        return genes

    def _arrow_head(self, x: float, y: float, strand: str) -> str:
        size = 6
        if strand == "+":
            points = [(x, y), (x - size, y - 4), (x - size, y + 4)]
        else:
            points = [(x, y), (x + size, y - 4), (x + size, y + 4)]
        point_text = " ".join(f"{px:.2f},{py:.2f}" for px, py in points)
        return f'<polygon points="{point_text}" fill="black" />'

    def _target_chromosome(self) -> str:
        return self._normalize_chromosome(str(self.args.chr))

    def _normalize_chromosome(self, value: str) -> str:
        value = str(value)
        if self.args.add_chr_prefix and not value.startswith(self.args.chr_prefix):
            return f"{self.args.chr_prefix}{value}"
        return value

    def _ensure_dirs(self) -> dict[str, Path]:
        dirs = {
            "base": self.output_dir,
            "plots": self.output_dir / "plots",
            "ld": self.output_dir / "ld",
            "logs": self.output_dir / "logs",
            "state": self.output_dir / "state",
        }
        for path in dirs.values():
            path.mkdir(parents=True, exist_ok=True)
        return dirs

    def _build_logger(self) -> logging.Logger:
        logger = logging.getLogger(f"gwas_step4_{id(self)}")
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
        description="GWAS step4 pipeline: local Manhattan plot and LD plot wrapper."
    )
    parser.add_argument("--gwas", required=True, help="GWAS result file")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--chr", required=True, help="Target chromosome")
    parser.add_argument("--start", required=True, type=int, help="Region start")
    parser.add_argument("--end", required=True, type=int, help="Region end")
    parser.add_argument("--chr-column", default="CHR", help="Chromosome column name")
    parser.add_argument("--pos-column", default="BP", help="Position column name")
    parser.add_argument("--p-column", default="P", help="P-value column name")
    parser.add_argument("--threshold-line", type=float, default=1e-5, help="Horizontal significance line")
    parser.add_argument("--gene-bed", help="Optional BED file for gene track")
    parser.add_argument("--vcf", help="Optional VCF file for LD analysis")
    parser.add_argument("--ldblockshow", help="Optional path to LDBlockShow executable")
    parser.add_argument("--add-chr-prefix", action="store_true", help="Add chromosome prefix before matching")
    parser.add_argument("--chr-prefix", default="Chr", help="Chromosome prefix when enabled")
    parser.add_argument("--force", action="store_true", help="Rerun even if outputs already exist")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    pipeline = Step4Pipeline(args)
    pipeline.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
