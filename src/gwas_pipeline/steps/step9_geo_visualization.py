from __future__ import annotations

import argparse
import csv
import json
import logging
import math
from collections import Counter, defaultdict
from pathlib import Path
import sys
from typing import Any


SAMPLE_TYPE_COLORS = {
    "Cultivar": "#7d3e87",
    "Landrace": "#fabb73",
    "Wild": "#8bccb7",
}

HAP_COLORS = {
    "Hap1": "#76c9b0",
    "Hap2": "#fd9972",
    "Hap3": "#99abd1",
    "Hap4": "#c7a0d8",
    "Hap5": "#f3de8a",
    "Hap6": "#8dd3c7",
    "Hap7": "#fb8072",
}


class Step9Pipeline:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.output_dir = Path(args.output_dir).resolve()
        self.dirs = self._ensure_dirs()
        self.logger = self._build_logger()
        self.state_path = self.dirs["state"] / "pipeline_state.json"
        self.state = self._load_state()

    def run(self) -> None:
        self.logger.info("Starting GWAS step9 geographic visualization workflow")
        samples = self.read_sample_points()
        hap_samples = self.read_haplotype_points()
        sample_site_table = self.aggregate_sample_sites(samples)
        hap_site_table = self.aggregate_haplotype_sites(hap_samples)
        regression_table = self.compute_latitude_regression(hap_samples)
        self.plot_sample_type_map(sample_site_table)
        self.plot_haplotype_map(hap_site_table)
        self.plot_latitude_regression(regression_table)
        self.logger.info("Step9 workflow finished successfully")

    def read_sample_points(self) -> list[dict[str, Any]]:
        if not self.args.samples:
            return []
        path = Path(self.args.samples).resolve()
        if not path.exists():
            raise FileNotFoundError(f"Missing sample file: {path}")
        points = []
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            for row in reader:
                if not row:
                    continue
                points.append(
                    {
                        "Type": row[0],
                        "Lat": float(row[1]),
                        "Lon": float(row[2]),
                        "City": row[3],
                    }
                )
        return points

    def read_haplotype_points(self) -> list[dict[str, Any]]:
        if not self.args.hap_samples:
            return []
        path = Path(self.args.hap_samples).resolve()
        if not path.exists():
            raise FileNotFoundError(f"Missing haplotype sample file: {path}")
        points = []
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            for row in reader:
                if not row:
                    continue
                points.append(
                    {
                        "Sample": row[0],
                        "Hap": row[1],
                        "Type": row[2],
                        "Lat": float(row[3]),
                        "Lon": float(row[4]),
                        "City": row[5],
                    }
                )
        return points

    def aggregate_sample_sites(self, samples: list[dict[str, Any]]) -> Path:
        out_path = self.output_dir / "sample_type_counts_by_site.tsv"
        if self._should_skip_global("aggregate_sample_sites", [out_path]):
            self.logger.info("Skipping sample-type site aggregation")
            return out_path

        counts: dict[tuple[float, float], Counter[str]] = defaultdict(Counter)
        for row in samples:
            counts[(row["Lon"], row["Lat"])][row["Type"]] += 1

        categories = self._ordered_sample_types(samples)
        rows = []
        for (lon, lat), counter in sorted(counts.items()):
            total = sum(counter.values())
            out = {
                "Lon": f"{lon:.8f}",
                "Lat": f"{lat:.8f}",
                "Total": str(total),
            }
            for cat in categories:
                out[cat] = str(counter.get(cat, 0))
            rows.append(out)

        self._write_tsv(out_path, ["Lon", "Lat", "Total"] + categories, rows)
        self._mark_global_done("aggregate_sample_sites")
        return out_path

    def aggregate_haplotype_sites(self, hap_samples: list[dict[str, Any]]) -> Path:
        out_path = self.output_dir / "haplotype_counts_by_site.tsv"
        if self._should_skip_global("aggregate_haplotype_sites", [out_path]):
            self.logger.info("Skipping haplotype site aggregation")
            return out_path

        counts: dict[tuple[float, float], Counter[str]] = defaultdict(Counter)
        for row in hap_samples:
            counts[(row["Lon"], row["Lat"])][row["Hap"]] += 1

        hap_names = self._ordered_haplotypes(hap_samples)
        rows = []
        for (lon, lat), counter in sorted(counts.items()):
            total = sum(counter.values())
            out = {
                "Lon": f"{lon:.8f}",
                "Lat": f"{lat:.8f}",
                "Total": str(total),
            }
            for hap in hap_names:
                out[hap] = str(counter.get(hap, 0))
            rows.append(out)

        self._write_tsv(out_path, ["Lon", "Lat", "Total"] + hap_names, rows)
        self._mark_global_done("aggregate_haplotype_sites")
        return out_path

    def compute_latitude_regression(self, hap_samples: list[dict[str, Any]]) -> Path:
        out_path = self.output_dir / "latitude_regression.tsv"
        if self._should_skip_global("compute_latitude_regression", [out_path]):
            self.logger.info("Skipping latitude regression summary")
            return out_path

        if not hap_samples:
            self._write_tsv(out_path, ["Hap", "Slope", "Intercept", "R2", "PValueApprox"], [])
            self._mark_global_done("compute_latitude_regression")
            return out_path

        bin_counts: dict[int, Counter[str]] = defaultdict(Counter)
        for row in hap_samples:
            abs_lat = abs(float(row["Lat"]))
            lat_bin = int(math.floor(abs_lat / 5) * 5 + 5)
            bin_counts[lat_bin][row["Hap"]] += 1

        hap_names = self._ordered_haplotypes(hap_samples)
        rows = []
        for hap in hap_names:
            x_vals = []
            y_vals = []
            for lat_bin in sorted(bin_counts):
                total = sum(bin_counts[lat_bin].values())
                proportion = bin_counts[lat_bin].get(hap, 0) / total if total else 0.0
                x_vals.append(proportion)
                y_vals.append(float(lat_bin))
            slope, intercept, r2 = self._linear_regression(x_vals, y_vals)
            p_approx = self._pvalue_label_from_r2(r2, len(x_vals))
            rows.append(
                {
                    "Hap": hap,
                    "Slope": f"{slope:.8g}",
                    "Intercept": f"{intercept:.8g}",
                    "R2": f"{r2:.6f}",
                    "PValueApprox": p_approx,
                }
            )

        self._write_tsv(out_path, ["Hap", "Slope", "Intercept", "R2", "PValueApprox"], rows)
        self._mark_global_done("compute_latitude_regression")
        return out_path

    def plot_sample_type_map(self, site_table: Path) -> Path:
        out_path = self.dirs["plots"] / "sample_type_map.svg"
        if self._should_skip_global("plot_sample_type_map", [out_path]):
            self.logger.info("Skipping sample-type map")
            return out_path

        rows = self._read_tsv_rows(site_table)
        categories = [key for key in rows[0].keys() if key not in {"Lon", "Lat", "Total"}] if rows else []
        svg = self._build_pie_map_svg(
            rows=rows,
            categories=categories,
            color_map={cat: SAMPLE_TYPE_COLORS.get(cat, "#999999") for cat in categories},
            title="Sample Type Geographic Distribution",
            legend_title="Type",
        )
        out_path.write_text(svg, encoding="utf-8")
        self._mark_global_done("plot_sample_type_map")
        return out_path

    def plot_haplotype_map(self, site_table: Path) -> Path:
        out_path = self.dirs["plots"] / "haplotype_map.svg"
        if self._should_skip_global("plot_haplotype_map", [out_path]):
            self.logger.info("Skipping haplotype map")
            return out_path

        rows = self._read_tsv_rows(site_table)
        categories = [key for key in rows[0].keys() if key not in {"Lon", "Lat", "Total"}] if rows else []
        svg = self._build_pie_map_svg(
            rows=rows,
            categories=categories,
            color_map={cat: HAP_COLORS.get(cat, "#999999") for cat in categories},
            title="Haplotype Geographic Distribution",
            legend_title="Haplotype",
        )
        out_path.write_text(svg, encoding="utf-8")
        self._mark_global_done("plot_haplotype_map")
        return out_path

    def plot_latitude_regression(self, regression_table: Path) -> Path:
        out_path = self.dirs["plots"] / "latitude_regression.svg"
        if self._should_skip_global("plot_latitude_regression", [out_path]):
            self.logger.info("Skipping latitude regression plot")
            return out_path

        if not self.args.hap_samples:
            out_path.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="800" height="500"></svg>', encoding="utf-8")
            self._mark_global_done("plot_latitude_regression")
            return out_path

        hap_samples = self.read_haplotype_points()
        bin_counts: dict[int, Counter[str]] = defaultdict(Counter)
        for row in hap_samples:
            abs_lat = abs(float(row["Lat"]))
            lat_bin = int(math.floor(abs_lat / 5) * 5 + 5)
            bin_counts[lat_bin][row["Hap"]] += 1

        hap_names = self._ordered_haplotypes(hap_samples)
        points_by_hap: dict[str, list[tuple[float, float]]] = {}
        for hap in hap_names:
            points = []
            for lat_bin in sorted(bin_counts):
                total = sum(bin_counts[lat_bin].values())
                proportion = bin_counts[lat_bin].get(hap, 0) / total if total else 0.0
                points.append((proportion, float(lat_bin)))
            points_by_hap[hap] = points

        stats = {row["Hap"]: row for row in self._read_tsv_rows(regression_table)}
        svg = self._build_regression_svg(points_by_hap, stats)
        out_path.write_text(svg, encoding="utf-8")
        self._mark_global_done("plot_latitude_regression")
        return out_path

    def _build_pie_map_svg(
        self,
        rows: list[dict[str, str]],
        categories: list[str],
        color_map: dict[str, str],
        title: str,
        legend_title: str,
    ) -> str:
        width = 1200
        height = 720
        left = 50
        top = 70
        map_width = 1000
        map_height = 580

        xlim = self.args.xlim if self.args.xlim else (-180.0, 180.0)
        ylim = self.args.ylim if self.args.ylim else (-60.0, 85.0)
        lon_min, lon_max = xlim
        lat_min, lat_max = ylim

        elements = [
            '<rect width="100%" height="100%" fill="white" />',
            f'<text x="{left}" y="34" font-size="28" font-weight="700">{self._escape(title)}</text>',
            f'<rect x="{left}" y="{top}" width="{map_width}" height="{map_height}" rx="18" fill="#dcdcdc" stroke="none" />',
        ]

        for lon in self._nice_ticks(lon_min, lon_max, 6):
            x = left + ((lon - lon_min) / max(1e-9, lon_max - lon_min)) * map_width
            elements.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + map_height}" stroke="#efefef" />')
            elements.append(f'<text x="{x:.2f}" y="{top + map_height + 22}" text-anchor="middle" font-size="11">{lon:g}</text>')

        for lat in self._nice_ticks(lat_min, lat_max, 6):
            y = top + map_height - ((lat - lat_min) / max(1e-9, lat_max - lat_min)) * map_height
            elements.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + map_width}" y2="{y:.2f}" stroke="#efefef" />')
            elements.append(f'<text x="{left - 10}" y="{y + 4:.2f}" text-anchor="end" font-size="11">{lat:g}</text>')

        max_total = max((int(row["Total"]) for row in rows), default=1)
        for row in rows:
            lon = float(row["Lon"])
            lat = float(row["Lat"])
            if not (lon_min <= lon <= lon_max and lat_min <= lat <= lat_max):
                continue
            total = int(row["Total"])
            cx = left + ((lon - lon_min) / max(1e-9, lon_max - lon_min)) * map_width
            cy = top + map_height - ((lat - lat_min) / max(1e-9, lat_max - lat_min)) * map_height
            radius = 5 + 16 * math.sqrt(total / max_total)
            elements.extend(self._pie_slices(cx, cy, radius, row, categories, color_map))

        legend_x = left + map_width + 35
        legend_y = top + 20
        elements.append(f'<text x="{legend_x}" y="{legend_y - 10}" font-size="15" font-weight="700">{self._escape(legend_title)}</text>')
        for idx, category in enumerate(categories):
            y = legend_y + idx * 24
            elements.append(f'<rect x="{legend_x}" y="{y}" width="16" height="16" rx="3" fill="{color_map.get(category, "#999")}" />')
            elements.append(f'<text x="{legend_x + 24}" y="{y + 13}" font-size="13">{self._escape(category)}</text>')

        return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">
{''.join(elements)}
</svg>
"""

    def _build_regression_svg(
        self,
        points_by_hap: dict[str, list[tuple[float, float]]],
        stats: dict[str, dict[str, str]],
    ) -> str:
        width = 900
        height = 900
        left = 90
        top = 70
        right = 50
        bottom = 80
        plot_width = width - left - right
        plot_height = height - top - bottom

        x_min = 0.0
        x_max = max((x for pts in points_by_hap.values() for x, _ in pts), default=1.0)
        x_max = max(1.0, x_max * 1.05)
        y_min = 0.0
        y_max = max((y for pts in points_by_hap.values() for _, y in pts), default=90.0)
        y_max = max(90.0, y_max)

        elements = [
            '<rect width="100%" height="100%" fill="white" />',
            f'<text x="{left}" y="34" font-size="28" font-weight="700">Haplotype Proportion vs Absolute Latitude</text>',
            f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#333" />',
            f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#333" />',
            f'<text x="{left + plot_width / 2:.2f}" y="{height - 18}" text-anchor="middle" font-size="16">Haplotype Proportion</text>',
            f'<text x="24" y="{top + plot_height / 2:.2f}" text-anchor="middle" font-size="16" transform="rotate(-90 24 {top + plot_height / 2:.2f})">Absolute Latitude</text>',
        ]

        for frac in [0.0, 0.25, 0.5, 0.75, 1.0]:
            x = left + frac * plot_width
            value = x_min + frac * (x_max - x_min)
            elements.append(f'<line x1="{x:.2f}" y1="{top + plot_height}" x2="{x:.2f}" y2="{top}" stroke="#f0f0f0" />')
            elements.append(f'<text x="{x:.2f}" y="{top + plot_height + 24:.2f}" text-anchor="middle" font-size="12">{value:.0%}</text>')
        for lat in range(0, int(y_max) + 1, 10):
            y = top + plot_height - ((lat - y_min) / max(1e-9, y_max - y_min)) * plot_height
            elements.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" stroke="#f0f0f0" />')
            elements.append(f'<text x="{left - 10}" y="{y + 4:.2f}" text-anchor="end" font-size="12">{lat}</text>')

        hap_names = list(points_by_hap.keys())
        color_map = {hap: HAP_COLORS.get(hap, "#999999") for hap in hap_names}
        for hap in hap_names:
            pts = points_by_hap[hap]
            coords = []
            for xval, yval in pts:
                x = left + ((xval - x_min) / max(1e-9, x_max - x_min)) * plot_width
                y = top + plot_height - ((yval - y_min) / max(1e-9, y_max - y_min)) * plot_height
                coords.append((x, y))
                elements.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4" fill="{color_map[hap]}" fill-opacity="0.85" />')

            slope = float(stats[hap]["Slope"])
            intercept = float(stats[hap]["Intercept"])
            x1_data = x_min
            x2_data = x_max
            y1_data = slope * x1_data + intercept
            y2_data = slope * x2_data + intercept
            x1 = left + ((x1_data - x_min) / max(1e-9, x_max - x_min)) * plot_width
            x2 = left + ((x2_data - x_min) / max(1e-9, x_max - x_min)) * plot_width
            y1 = top + plot_height - ((y1_data - y_min) / max(1e-9, y_max - y_min)) * plot_height
            y2 = top + plot_height - ((y2_data - y_min) / max(1e-9, y_max - y_min)) * plot_height
            elements.append(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color_map[hap]}" stroke-width="2.2" />')

        label_x = left + plot_width - 10
        label_y = top + 20
        for idx, hap in enumerate(hap_names):
            stat = stats[hap]
            y = label_y + idx * 24
            label = f"{hap}: R2={float(stat['R2']):.3f}, p={stat['PValueApprox']}"
            elements.append(f'<text x="{label_x}" y="{y}" text-anchor="end" font-size="13" fill="{color_map[hap]}">{self._escape(label)}</text>')

        return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">
{''.join(elements)}
</svg>
"""

    def _pie_slices(
        self,
        cx: float,
        cy: float,
        radius: float,
        row: dict[str, str],
        categories: list[str],
        color_map: dict[str, str],
    ) -> list[str]:
        total = sum(int(row.get(cat, 0)) for cat in categories)
        if total == 0:
            return []
        parts = []
        angle = -math.pi / 2
        for cat in categories:
            count = int(row.get(cat, 0))
            if count <= 0:
                continue
            delta = 2 * math.pi * (count / total)
            parts.append(self._pie_slice(cx, cy, radius, angle, angle + delta, color_map.get(cat, "#999999")))
            angle += delta
        parts.append(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{radius:.2f}" fill="none" stroke="#ffffff" stroke-width="0.8" />')
        return parts

    def _pie_slice(self, cx: float, cy: float, radius: float, start_angle: float, end_angle: float, color: str) -> str:
        x1 = cx + radius * math.cos(start_angle)
        y1 = cy + radius * math.sin(start_angle)
        x2 = cx + radius * math.cos(end_angle)
        y2 = cy + radius * math.sin(end_angle)
        large_arc = 1 if end_angle - start_angle > math.pi else 0
        return (
            f'<path d="M {cx:.2f} {cy:.2f} L {x1:.2f} {y1:.2f} '
            f'A {radius:.2f} {radius:.2f} 0 {large_arc} 1 {x2:.2f} {y2:.2f} Z" '
            f'fill="{color}" fill-opacity="0.95" />'
        )

    def _ordered_sample_types(self, samples: list[dict[str, Any]]) -> list[str]:
        preferred = ["Cultivar", "Landrace", "Wild"]
        seen = {row["Type"] for row in samples}
        ordered = [item for item in preferred if item in seen]
        ordered.extend(sorted(item for item in seen if item not in ordered))
        return ordered

    def _ordered_haplotypes(self, hap_samples: list[dict[str, Any]]) -> list[str]:
        seen = sorted({row["Hap"] for row in hap_samples})
        return seen

    def _linear_regression(self, x_vals: list[float], y_vals: list[float]) -> tuple[float, float, float]:
        n = len(x_vals)
        if n == 0:
            return 0.0, 0.0, 0.0
        mean_x = sum(x_vals) / n
        mean_y = sum(y_vals) / n
        ss_xx = sum((x - mean_x) ** 2 for x in x_vals)
        ss_xy = sum((x - mean_x) * (y - mean_y) for x, y in zip(x_vals, y_vals))
        slope = ss_xy / ss_xx if ss_xx else 0.0
        intercept = mean_y - slope * mean_x
        y_hat = [slope * x + intercept for x in x_vals]
        ss_tot = sum((y - mean_y) ** 2 for y in y_vals)
        ss_res = sum((y - yh) ** 2 for y, yh in zip(y_vals, y_hat))
        r2 = 1.0 - ss_res / ss_tot if ss_tot else 0.0
        return slope, intercept, max(0.0, r2)

    def _pvalue_label_from_r2(self, r2: float, n: int) -> str:
        if n < 3:
            return "NA"
        if r2 >= 0.60:
            return "<0.001"
        if r2 >= 0.35:
            return "0.001-0.01"
        if r2 >= 0.15:
            return "0.01-0.05"
        return ">0.05"

    def _nice_ticks(self, min_val: float, max_val: float, n: int) -> list[float]:
        if n <= 1:
            return [min_val, max_val]
        step = (max_val - min_val) / (n - 1)
        return [round(min_val + i * step, 2) for i in range(n)]

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
        logger = logging.getLogger(f"gwas_step9_{id(self)}")
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
        description="GWAS step9 pipeline: sample-type map, haplotype map, and latitude regression."
    )
    parser.add_argument("--samples", help="Headerless sample-type CSV")
    parser.add_argument("--hap-samples", help="Headerless haplotype sample CSV")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--xlim", nargs=2, type=float, help="Optional longitude range, for example -10 50")
    parser.add_argument("--ylim", nargs=2, type=float, help="Optional latitude range, for example 20 70")
    parser.add_argument("--force", action="store_true", help="Rerun even if outputs already exist")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.samples and not args.hap_samples:
        raise SystemExit("Provide at least one of --samples or --hap-samples")
    pipeline = Step9Pipeline(args)
    pipeline.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
