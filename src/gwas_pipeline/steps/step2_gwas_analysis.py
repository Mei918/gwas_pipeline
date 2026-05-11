from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path
import shlex
import subprocess
import sys
from typing import Any, Optional


class Step2Pipeline:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.output_dir = Path(args.output_dir).resolve()
        self.dirs = self._ensure_dirs()
        self.logger = self._build_logger()
        self.state_path = self.dirs["state"] / "pipeline_state.json"
        self.state = self._load_state()

    def run(self) -> None:
        if self.args.mode == "gcta":
            self.run_gcta_mode()
        elif self.args.mode == "rmvp":
            self.run_rmvp_mode()
        else:
            raise ValueError(f"Unsupported mode: {self.args.mode}")

    def run_gcta_mode(self) -> None:
        self.logger.info("Starting GCTA MLM workflow")
        self._check_required_args(["bfile", "phenotype"])

        phenotype_in = Path(self.args.phenotype).resolve()
        bfile = Path(self.args.bfile).resolve()

        prepared_pheno = self.prepare_gcta_phenotype(phenotype_in)
        prepared_fam_prefix = self.copy_plink_with_updated_fam(bfile, prepared_pheno)
        grm_prefix = self.make_grm(prepared_fam_prefix)
        pca_prefix = self.make_pca(prepared_fam_prefix)
        result_file = self.run_gcta_mlm(prepared_fam_prefix, prepared_pheno, grm_prefix, pca_prefix)
        self.make_plots_from_mlma(result_file)
        self.logger.info("GCTA MLM workflow finished successfully")

    def run_rmvp_mode(self) -> None:
        self.logger.info("Starting rMVP workflow")
        self._check_required_args(["vcf", "phenotype"])

        vcf = Path(self.args.vcf).resolve()
        phenotype_in = Path(self.args.phenotype).resolve()

        prepared_pheno = self.prepare_rmvp_phenotype(phenotype_in)
        plink_prefix = self.make_rmvp_support_files(vcf)
        r_script = self.write_rmvp_script(vcf, prepared_pheno, plink_prefix)
        self.run_r_script(r_script)
        self.logger.info("rMVP workflow finished successfully")

    def prepare_gcta_phenotype(self, phenotype_in: Path) -> Path:
        out_path = self.output_dir / "phenotype.txt"
        if self._should_skip_global("prepare_gcta_phenotype", [out_path]):
            self.logger.info("Skipping GCTA phenotype preparation")
            return out_path

        rows = []
        with phenotype_in.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                fields = line.split()
                if len(fields) < 2:
                    raise ValueError(f"Invalid phenotype line: {line}")
                sample_id = fields[0]
                phenotype = fields[1]
                rows.append(f"{sample_id}\t{sample_id}\t{phenotype}")

        out_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
        self._mark_global_done("prepare_gcta_phenotype")
        return out_path

    def prepare_rmvp_phenotype(self, phenotype_in: Path) -> Path:
        out_path = self.output_dir / "phenotype_rmvp.txt"
        if self._should_skip_global("prepare_rmvp_phenotype", [out_path]):
            self.logger.info("Skipping rMVP phenotype preparation")
            return out_path

        rows = []
        with phenotype_in.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                fields = line.split()
                if len(fields) < 2:
                    raise ValueError(f"Invalid phenotype line: {line}")
                rows.append(f"{fields[0]}\t{fields[1]}")

        out_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
        self._mark_global_done("prepare_rmvp_phenotype")
        return out_path

    def copy_plink_with_updated_fam(self, source_prefix: Path, phenotype_txt: Path) -> Path:
        target_prefix = self.output_dir / "cohort_working"
        fam_out = Path(f"{target_prefix}.fam")
        bed_out = Path(f"{target_prefix}.bed")
        bim_out = Path(f"{target_prefix}.bim")

        if self._should_skip_global("prepare_working_bfile", [fam_out, bed_out, bim_out]):
            self.logger.info("Skipping working PLINK copy")
            return target_prefix

        for ext in ["bed", "bim", "fam"]:
            src = Path(f"{source_prefix}.{ext}")
            dst = Path(f"{target_prefix}.{ext}")
            if not src.exists():
                raise FileNotFoundError(f"Missing PLINK input file: {src}")
            dst.write_bytes(src.read_bytes())

        phenotype_by_id: dict[str, str] = {}
        with phenotype_txt.open("r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.strip().split()
                if len(fields) >= 3:
                    phenotype_by_id[fields[1]] = fields[2]

        updated_lines = []
        with fam_out.open("r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split()
                if len(fields) < 6:
                    raise ValueError(f"Invalid FAM line: {line}")
                iid = fields[1]
                if iid in phenotype_by_id:
                    fields[5] = phenotype_by_id[iid]
                updated_lines.append("\t".join(fields))

        fam_out.write_text("\n".join(updated_lines) + "\n", encoding="utf-8")
        self._mark_global_done("prepare_working_bfile")
        return target_prefix

    def make_grm(self, bfile_prefix: Path) -> Path:
        out_prefix = self.output_dir / "kinship_matrix"
        required = [Path(f"{out_prefix}.grm.bin"), Path(f"{out_prefix}.grm.N.bin"), Path(f"{out_prefix}.grm.id")]
        if self._should_skip_global("make_grm", required):
            self.logger.info("Skipping GRM construction")
            return out_prefix

        cmd = f"gcta64 --bfile {self._q(bfile_prefix)} --make-grm --out {self._q(out_prefix)}"
        self._run_cmd(cmd)
        self._mark_global_done("make_grm")
        return out_prefix

    def make_pca(self, bfile_prefix: Path) -> Path:
        out_prefix = self.output_dir / "gwas_pca"
        required = [Path(f"{out_prefix}.eigenvec"), Path(f"{out_prefix}.eigenval")]
        if self._should_skip_global("make_pca", required):
            self.logger.info("Skipping PCA")
            return out_prefix

        cmd = (
            f"plink --bfile {self._q(bfile_prefix)} "
            f"--pca {self.args.pca_components} "
            f"--out {self._q(out_prefix)}"
        )
        self._run_cmd(cmd)
        self._mark_global_done("make_pca")
        return out_prefix

    def run_gcta_mlm(
        self,
        bfile_prefix: Path,
        phenotype_txt: Path,
        grm_prefix: Path,
        pca_prefix: Path,
    ) -> Path:
        out_prefix = self.output_dir / "gwas_mlm_results"
        mlma_file = Path(f"{out_prefix}.mlma")
        if self._should_skip_global("run_gcta_mlm", [mlma_file]):
            self.logger.info("Skipping GCTA MLM analysis")
            return mlma_file

        cmd = (
            f"gcta64 --bfile {self._q(bfile_prefix)} "
            f"--mlma "
            f"--pheno {self._q(phenotype_txt)} "
            f"--qcovar {self._q(Path(f'{pca_prefix}.eigenvec'))} "
            f"--grm {self._q(grm_prefix)} "
            f"--out {self._q(out_prefix)}"
        )
        self._run_cmd(cmd)
        self._mark_global_done("run_gcta_mlm")
        return mlma_file

    def make_rmvp_support_files(self, vcf: Path) -> Path:
        out_prefix = self.output_dir / "gwas_output"
        required = [Path(f"{out_prefix}.eigenvec")]
        if self._should_skip_global("make_rmvp_support_files", required):
            self.logger.info("Skipping rMVP support file generation")
            return out_prefix

        king_cmd = (
            f"plink --vcf {self._q(vcf)} "
            f"--make-king-table "
            f"--out {self._q(out_prefix)}"
        )
        pca_cmd = (
            f"plink --vcf {self._q(vcf)} "
            f"--pca 10 "
            f"--out {self._q(out_prefix)}"
        )
        self._run_cmd(king_cmd)
        self._run_cmd(pca_cmd)
        self._mark_global_done("make_rmvp_support_files")
        return out_prefix

    def write_rmvp_script(self, vcf: Path, phenotype_txt: Path, out_prefix: Path) -> Path:
        script_path = self.output_dir / "run_rmvp.R"
        if self._should_skip_global("write_rmvp_script", [script_path]):
            self.logger.info("Skipping rMVP script generation")
            return script_path

        content = f"""library(rMVP)
library(bigmemory)

MVP.Data(
  fileVCF = "{self._r_quote(vcf)}",
  filePhe = "{self._r_quote(phenotype_txt)}",
  sep.phe = "\\t",
  fileKin = TRUE,
  filePC = TRUE,
  out = "{self._r_quote(out_prefix)}"
)

genotype <- attach.big.matrix("{self._r_quote(Path(str(out_prefix) + '.geno.desc'))}")
phenotype <- read.table("{self._r_quote(Path(str(out_prefix) + '.phe'))}", header = TRUE)
map <- read.table("{self._r_quote(Path(str(out_prefix) + '.geno.map'))}", header = TRUE)
Kinship <- attach.big.matrix("{self._r_quote(Path(str(out_prefix) + '.kin.desc'))}")
Covariates_PC <- bigmemory::as.matrix(attach.big.matrix("{self._r_quote(Path(str(out_prefix) + '.pc.desc'))}"))

result <- MVP(
  phe = phenotype,
  geno = genotype,
  map = map,
  K = Kinship,
  CV.MLM = Covariates_PC,
  nPC.MLM = 0,
  vc.method = "BRENT",
  maxLoop = 10,
  threshold = 0.05,
  method = "MLM",
  file.output = c("pmap", "pmap.signal", "plot", "log")
)

write.table(result$pmap, file = "{self._r_quote(self.output_dir / 'rmvp_pmap.tsv')}", sep = "\\t", row.names = FALSE, quote = FALSE)
"""
        script_path.write_text(content, encoding="utf-8")
        self._mark_global_done("write_rmvp_script")
        return script_path

    def run_r_script(self, script_path: Path) -> None:
        marker = self.output_dir / "rmvp_pmap.tsv"
        if self._should_skip_global("run_r_script", [marker]):
            self.logger.info("Skipping rMVP execution")
            return

        cmd = f"Rscript {self._q(script_path)}"
        self._run_cmd(cmd)
        self._mark_global_done("run_r_script")

    def make_plots_from_mlma(self, mlma_file: Path) -> None:
        manhattan_svg = self.dirs["plots"] / "manhattan.svg"
        qq_svg = self.dirs["plots"] / "qq.svg"
        if self._should_skip_global("make_plots", [manhattan_svg, qq_svg]):
            self.logger.info("Skipping plot generation")
            return

        rows = self._read_mlma_rows(mlma_file)
        if not rows:
            raise RuntimeError(f"No GWAS rows found in {mlma_file}")

        self._write_manhattan_svg(manhattan_svg, rows)
        self._write_qq_svg(qq_svg, rows)
        self._mark_global_done("make_plots")

    def _read_mlma_rows(self, mlma_file: Path) -> list[dict[str, Any]]:
        rows = []
        with mlma_file.open("r", encoding="utf-8") as handle:
            header = handle.readline().strip().split()
            for line in handle:
                fields = line.strip().split()
                if len(fields) != len(header):
                    continue
                row = dict(zip(header, fields))
                p_value = float(row["p"])
                if p_value <= 0:
                    p_value = 1e-300
                rows.append(
                    {
                        "Chr": row["Chr"],
                        "SNP": row["SNP"],
                        "bp": int(row["bp"]),
                        "p": p_value,
                    }
                )
        return rows

    def _write_manhattan_svg(self, path: Path, rows: list[dict[str, Any]]) -> None:
        width = 1200
        height = 700
        left = 80
        right = 30
        top = 70
        bottom = 70
        plot_width = width - left - right
        plot_height = height - top - bottom

        ordered = sorted(rows, key=lambda row: (self._chr_key(row["Chr"]), row["bp"]))
        xs = list(range(len(ordered)))
        ys = [-math.log10(row["p"]) for row in ordered]
        max_y = max(ys) * 1.08 if ys else 1.0
        if max_y < 1.0:
            max_y = 1.0

        colors = ["#7a983e", "#b8cc8e"]
        circles = []
        chr_centers: dict[str, list[int]] = {}
        for idx, row in enumerate(ordered):
            chr_name = str(row["Chr"])
            chr_centers.setdefault(chr_name, []).append(idx)
            x = left + (idx / max(1, len(ordered) - 1)) * plot_width
            y = top + plot_height - (ys[idx] / max_y) * plot_height
            color = colors[self._chr_key(chr_name)[0] % 2]
            circles.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.8" fill="{color}" />')

        axis = [
            f'<line x1="{left}" y1="{top + plot_height}" x2="{width - right}" y2="{top + plot_height}" stroke="#333" />',
            f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#333" />',
            f'<text x="{width / 2:.0f}" y="32" text-anchor="middle" font-size="26">Manhattan Plot</text>',
            f'<text x="{width / 2:.0f}" y="{height - 18}" text-anchor="middle" font-size="16">Chromosome</text>',
            f'<text x="24" y="{height / 2:.0f}" text-anchor="middle" font-size="16" transform="rotate(-90 24 {height / 2:.0f})">-log10(p)</text>',
        ]
        for frac in [0.25, 0.5, 0.75, 1.0]:
            y = top + plot_height - frac * plot_height
            value = frac * max_y
            axis.append(f'<line x1="{left}" y1="{y:.2f}" x2="{width - right}" y2="{y:.2f}" stroke="#ececec" />')
            axis.append(f'<text x="{left - 10}" y="{y + 5:.2f}" text-anchor="end" font-size="12">{value:.1f}</text>')
        for chr_name, indexes in chr_centers.items():
            center_idx = (indexes[0] + indexes[-1]) / 2
            x = left + (center_idx / max(1, len(ordered) - 1)) * plot_width
            axis.append(f'<text x="{x:.2f}" y="{height - 40}" text-anchor="middle" font-size="12">{chr_name}</text>')

        svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">
<rect width="100%" height="100%" fill="white" />
{''.join(axis)}
{''.join(circles)}
</svg>
"""
        path.write_text(svg, encoding="utf-8")

    def _write_qq_svg(self, path: Path, rows: list[dict[str, Any]]) -> None:
        width = 700
        height = 800
        left = 80
        right = 30
        top = 70
        bottom = 70
        plot_width = width - left - right
        plot_height = height - top - bottom

        p_values = sorted(max(row["p"], 1e-300) for row in rows)
        expected = [-math.log10((i + 1) / (len(p_values) + 1)) for i in range(len(p_values))]
        observed = [-math.log10(p) for p in p_values]
        max_axis = max(max(expected), max(observed)) * 1.05
        if max_axis < 1.0:
            max_axis = 1.0

        points = []
        for exp, obs in zip(expected, observed):
            x = left + (exp / max_axis) * plot_width
            y = top + plot_height - (obs / max_axis) * plot_height
            points.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.5" fill="#5c7cfa" />')

        diagonal = (
            f'<line x1="{left}" y1="{top + plot_height}" '
            f'x2="{left + plot_width}" y2="{top}" stroke="#666" stroke-dasharray="6,4" />'
        )
        axis = [
            f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#333" />',
            f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#333" />',
            f'<text x="{width / 2:.0f}" y="32" text-anchor="middle" font-size="26">QQ Plot</text>',
            f'<text x="{width / 2:.0f}" y="{height - 18}" text-anchor="middle" font-size="16">Expected -log10(p)</text>',
            f'<text x="24" y="{height / 2:.0f}" text-anchor="middle" font-size="16" transform="rotate(-90 24 {height / 2:.0f})">Observed -log10(p)</text>',
        ]

        svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">
<rect width="100%" height="100%" fill="white" />
{''.join(axis)}
{diagonal}
{''.join(points)}
</svg>
"""
        path.write_text(svg, encoding="utf-8")

    def _ensure_dirs(self) -> dict[str, Path]:
        dirs = {
            "base": self.output_dir,
            "logs": self.output_dir / "logs",
            "state": self.output_dir / "state",
            "plots": self.output_dir / "plots",
        }
        for path in dirs.values():
            path.mkdir(parents=True, exist_ok=True)
        return dirs

    def _build_logger(self) -> logging.Logger:
        logger = logging.getLogger(f"gwas_step2_{id(self)}")
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

    def _check_required_args(self, names: list[str]) -> None:
        for name in names:
            if getattr(self.args, name) is None:
                raise ValueError(f"--{name.replace('_', '-')} is required in {self.args.mode} mode")

    @staticmethod
    def _q(path: Path) -> str:
        return shlex.quote(str(path))

    @staticmethod
    def _r_quote(path: Path) -> str:
        return str(path).replace("\\", "/")

    @staticmethod
    def _chr_key(value: str) -> tuple[int, str]:
        try:
            return (int(value), value)
        except ValueError:
            return (10**9, value)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="GWAS step2 pipeline: GCTA MLM or rMVP workflow with plotting."
    )
    parser.add_argument("--mode", choices=["gcta", "rmvp"], required=True, help="Workflow mode")
    parser.add_argument("--bfile", help="PLINK prefix for GCTA mode")
    parser.add_argument("--vcf", help="VCF path for rMVP mode")
    parser.add_argument("--phenotype", required=True, help="Phenotype file path")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--pca-components", type=int, default=5, help="Number of PCA components for GCTA mode")
    parser.add_argument("--force", action="store_true", help="Rerun even if outputs already exist")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    pipeline = Step2Pipeline(args)
    pipeline.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
