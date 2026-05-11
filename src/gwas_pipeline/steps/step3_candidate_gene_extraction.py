from __future__ import annotations

import argparse
import json
import logging
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Any, Optional


@dataclass
class GeneInterval:
    chrom: str
    start0: int
    end1: int
    gene_id: str
    strand: str


class Step3Pipeline:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.output_dir = Path(args.output_dir).resolve()
        self.dirs = self._ensure_dirs()
        self.logger = self._build_logger()
        self.state_path = self.dirs["state"] / "pipeline_state.json"
        self.state = self._load_state()

    def run(self) -> None:
        self.logger.info("Starting GWAS step3 candidate gene extraction")
        significant_hits = self.extract_significant_hits()
        gene_bed = self.build_gene_promoter_bed()
        annotated_hits = self.annotate_hits_with_genes(significant_hits, gene_bed)
        self.write_candidate_gene_summary(annotated_hits)
        self.logger.info("Step3 workflow finished successfully")

    def extract_significant_hits(self) -> Path:
        out_path = self.output_dir / "significant_hits.tsv"
        if self._should_skip_global("extract_significant_hits", [out_path]):
            self.logger.info("Skipping significant hit extraction")
            return out_path

        gwas_path = Path(self.args.gwas).resolve()
        with gwas_path.open("r", encoding="utf-8") as handle:
            header_line = handle.readline().strip()
            if not header_line:
                raise ValueError("GWAS file is empty")
            headers = header_line.split()

            p_column = self.args.p_column
            if p_column not in headers:
                raise ValueError(f"P-value column '{p_column}' not found in GWAS file")
            p_index = headers.index(p_column)

            out_lines = ["\t".join(headers)]
            kept = 0
            for line in handle:
                fields = line.strip().split()
                if len(fields) != len(headers):
                    continue
                p_value = float(fields[p_index])
                if p_value < self.args.p_threshold:
                    out_lines.append("\t".join(fields))
                    kept += 1

        out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
        self.logger.info("Retained %d significant hits at p < %s", kept, self.args.p_threshold)
        self._mark_global_done("extract_significant_hits")
        return out_path

    def build_gene_promoter_bed(self) -> Path:
        out_path = self.output_dir / "genes_with_promoter.bed"
        if self._should_skip_global("build_gene_promoter_bed", [out_path]):
            self.logger.info("Skipping gene+promoter BED generation")
            return out_path

        gff_path = Path(self.args.gff).resolve()
        intervals = []
        with gff_path.open("r", encoding="utf-8") as handle:
            for raw in handle:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 9:
                    continue
                chrom, _, feature_type, start, end, _, strand, _, attributes = fields
                if feature_type != "gene":
                    continue

                start1 = int(start)
                end1 = int(end)
                gene_id = self._extract_gene_id(attributes)
                if not gene_id:
                    continue

                if strand == "+":
                    ext_start1 = max(1, start1 - self.args.promoter_length)
                    ext_end1 = end1
                elif strand == "-":
                    ext_start1 = start1
                    ext_end1 = end1 + self.args.promoter_length
                else:
                    ext_start1 = start1
                    ext_end1 = end1

                interval = GeneInterval(
                    chrom=chrom,
                    start0=ext_start1 - 1,
                    end1=ext_end1,
                    gene_id=gene_id,
                    strand=strand,
                )
                intervals.append(interval)

        with out_path.open("w", encoding="utf-8") as handle:
            for interval in intervals:
                handle.write(
                    f"{interval.chrom}\t{interval.start0}\t{interval.end1}\t"
                    f"{interval.gene_id}\t.\t{interval.strand}\n"
                )

        self.logger.info("Wrote %d gene+promoter intervals", len(intervals))
        self._mark_global_done("build_gene_promoter_bed")
        return out_path

    def annotate_hits_with_genes(self, significant_hits: Path, gene_bed: Path) -> Path:
        out_path = self.output_dir / "significant_hits.with_gene.tsv"
        if self._should_skip_global("annotate_hits_with_genes", [out_path]):
            self.logger.info("Skipping significant hit annotation")
            return out_path

        intervals = self._load_gene_intervals(gene_bed)

        with significant_hits.open("r", encoding="utf-8") as handle:
            header_line = handle.readline().strip()
            if not header_line:
                raise ValueError("Significant hit file is empty")
            headers = header_line.split("\t")

            if self.args.chr_column not in headers:
                raise ValueError(f"Chromosome column '{self.args.chr_column}' not found")
            if self.args.pos_column not in headers:
                raise ValueError(f"Position column '{self.args.pos_column}' not found")

            chr_idx = headers.index(self.args.chr_column)
            pos_idx = headers.index(self.args.pos_column)

            out_lines = ["\t".join(headers + ["Gene"])]
            annotated_count = 0

            for raw in handle:
                fields = raw.strip().split("\t")
                if len(fields) != len(headers):
                    fields = raw.strip().split()
                if len(fields) != len(headers):
                    continue

                snp_chr = self._normalize_chromosome(fields[chr_idx])
                snp_pos = int(fields[pos_idx])
                hits = []
                for interval in intervals.get(snp_chr, []):
                    if interval.start0 <= snp_pos - 1 < interval.end1:
                        hits.append(interval.gene_id)

                if hits:
                    annotated_count += 1
                    out_lines.append("\t".join(fields + [",".join(hits)]))

        out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
        self.logger.info("Annotated %d significant hits to genes/promoters", annotated_count)
        self._mark_global_done("annotate_hits_with_genes")
        return out_path

    def write_candidate_gene_summary(self, annotated_hits_path: Path) -> Path:
        out_path = self.output_dir / "candidate_gene_summary.tsv"
        if self._should_skip_global("write_candidate_gene_summary", [out_path]):
            self.logger.info("Skipping candidate gene summary")
            return out_path

        gene_counts: dict[str, int] = {}
        with annotated_hits_path.open("r", encoding="utf-8") as handle:
            header = handle.readline()
            for raw in handle:
                fields = raw.strip().split("\t")
                if not fields:
                    continue
                genes = fields[-1].split(",")
                for gene in genes:
                    gene_counts[gene] = gene_counts.get(gene, 0) + 1

        rows = ["Gene\tHitCount"]
        for gene, count in sorted(gene_counts.items(), key=lambda item: (-item[1], item[0])):
            rows.append(f"{gene}\t{count}")
        out_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

        self.logger.info("Summarized %d candidate genes", len(gene_counts))
        self._mark_global_done("write_candidate_gene_summary")
        return out_path

    def _load_gene_intervals(self, gene_bed: Path) -> dict[str, list[GeneInterval]]:
        intervals: dict[str, list[GeneInterval]] = {}
        with gene_bed.open("r", encoding="utf-8") as handle:
            for raw in handle:
                fields = raw.strip().split("\t")
                if len(fields) < 6:
                    continue
                chrom = self._normalize_chromosome(fields[0])
                interval = GeneInterval(
                    chrom=chrom,
                    start0=int(fields[1]),
                    end1=int(fields[2]),
                    gene_id=fields[3],
                    strand=fields[5],
                )
                intervals.setdefault(chrom, []).append(interval)
        return intervals

    def _extract_gene_id(self, attributes: str) -> str:
        for item in attributes.split(";"):
            item = item.strip()
            if item.startswith("ID="):
                return item[3:]
        return ""

    def _normalize_chromosome(self, value: str) -> str:
        if self.args.add_chr_prefix and not value.startswith(self.args.chr_prefix):
            return f"{self.args.chr_prefix}{value}"
        return value

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
        logger = logging.getLogger(f"gwas_step3_{id(self)}")
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
        description="GWAS step3 pipeline: significant hit extraction and candidate gene annotation."
    )
    parser.add_argument("--gwas", required=True, help="GWAS result file")
    parser.add_argument("--gff", required=True, help="Gene annotation GFF")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--p-threshold", type=float, default=1e-5, help="P-value threshold")
    parser.add_argument("--p-column", default="P", help="P-value column name in GWAS file")
    parser.add_argument("--chr-column", default="CHR", help="Chromosome column name in GWAS file")
    parser.add_argument("--pos-column", default="BP", help="Position column name in GWAS file")
    parser.add_argument("--promoter-length", type=int, default=2000, help="Promoter extension length")
    parser.add_argument(
        "--add-chr-prefix",
        action="store_true",
        help="Add chromosome prefix to GWAS chromosome values before matching",
    )
    parser.add_argument(
        "--chr-prefix",
        default="Chr",
        help="Chromosome prefix used when --add-chr-prefix is enabled",
    )
    parser.add_argument("--force", action="store_true", help="Rerun even if outputs already exist")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    pipeline = Step3Pipeline(args)
    pipeline.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
