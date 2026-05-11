from __future__ import annotations

import argparse
import csv
import json
import logging
from dataclasses import dataclass
from pathlib import Path
import re
import shlex
import subprocess
import sys
from typing import Any, Optional


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMPLEMENT = str.maketrans({"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"})
VARIANT_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<pos>\d+)-(?P<ref>[A-Z]+)/(?P<alt>[A-Z]+)$")


@dataclass
class CdsGene:
    gene_id: str
    chrom: str
    strand: str
    cds_regions: list[tuple[int, int]]


class Step5Pipeline:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.output_dir = Path(args.output_dir).resolve()
        self.dirs = self._ensure_dirs()
        self.logger = self._build_logger()
        self.state_path = self.dirs["state"] / "pipeline_state.json"
        self.state = self._load_state()

    def run(self) -> None:
        self.logger.info("Starting GWAS step5 amino-acid mutation workflow")
        region_vcf = self.prepare_region_vcf()
        csv_path = self.vcf_to_csv(region_vcf)
        gff_path = Path(self.args.gff).resolve()
        cds_fasta = Path(self.args.cds_fasta).resolve()
        mut_csv = self.annotate_variant_effects(csv_path, gff_path, cds_fasta)
        self.write_effect_summary(mut_csv)
        self.logger.info("Step5 workflow finished successfully")

    def prepare_region_vcf(self) -> Path:
        if self.args.region_vcf:
            path = Path(self.args.region_vcf).resolve()
            if not path.exists():
                raise FileNotFoundError(f"Missing region VCF: {path}")
            return path

        if not self.args.vcf or not self.args.region:
            raise ValueError("Provide either --region-vcf or both --vcf and --region")

        out_path = self.output_dir / "gene_region.vcf"
        if self._should_skip_global("extract_region_vcf", [out_path]):
            self.logger.info("Skipping regional VCF extraction")
            return out_path

        vcf = Path(self.args.vcf).resolve()
        cmd = (
            f"bcftools view {self._q(vcf)} "
            f"-r {shlex.quote(self.args.region)} "
            f"-Ov -o {self._q(out_path)}"
        )
        self._run_cmd(cmd)
        self._mark_global_done("extract_region_vcf")
        return out_path

    def vcf_to_csv(self, region_vcf: Path) -> Path:
        out_path = self.output_dir / "gene_region.csv"
        if self._should_skip_global("vcf_to_csv", [out_path]):
            self.logger.info("Skipping VCF to CSV conversion")
            return out_path

        with region_vcf.open("r", encoding="utf-8") as handle:
            lines = [line.rstrip("\n") for line in handle if not line.startswith("##")]

        if not lines:
            raise ValueError("Region VCF has no readable lines")

        header = lines[0].split()
        samples = header[9:]
        variants = []
        haplotypes = {sample: [] for sample in samples}

        for line in lines[1:]:
            fields = line.split()
            if len(fields) < 10:
                continue
            chrom, pos, _, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            variants.append(f"{chrom}:{pos}-{ref}/{alt}")

            for idx, sample in enumerate(samples):
                gt_full = fields[9 + idx]
                gt = gt_full.split(":", 1)[0]
                if "." in gt:
                    haplotypes[sample].append("./.")
                    continue
                separator = "/" if "/" in gt else "|"
                alleles = gt.split(separator)
                bases = []
                for allele in alleles:
                    if allele == "0":
                        bases.append(ref)
                    elif allele == "1":
                        bases.append(alt)
                    else:
                        bases.append("N")
                haplotypes[sample].append("/".join(bases))

        with out_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["Sample"] + variants)
            for sample in samples:
                writer.writerow([sample] + haplotypes[sample])

        self._mark_global_done("vcf_to_csv")
        return out_path

    def annotate_variant_effects(self, snp_csv: Path, gff_file: Path, cds_fasta: Path) -> Path:
        out_path = self.output_dir / "gene_region_mut.csv"
        summary_path = self.output_dir / "variant_effect_summary.tsv"
        if self._should_skip_global("annotate_variant_effects", [out_path, summary_path]):
            self.logger.info("Skipping amino-acid effect annotation")
            return out_path

        cds_map = self.parse_gff(gff_file)
        cds_sequences = self.parse_fasta(cds_fasta)
        header, rows = self.read_csv_matrix(snp_csv)

        variant_columns = header[1:]
        results_by_column: dict[str, dict[str, str]] = {}
        summary_rows = [["Original_Column", "Gene_ID", "CDS_Position", "Original_Codon", "Mutated_Codon", "Effect"]]

        for col in variant_columns:
            effect = self.annotate_single_variant(col, cds_map, cds_sequences)
            if effect is not None:
                results_by_column[col] = effect
                summary_rows.append([
                    col,
                    effect["Gene_ID"],
                    str(effect["CDS_Position"]),
                    effect["Original_Codon"],
                    effect["Mutated_Codon"],
                    effect["Effect"],
                ])

        new_header = ["Sample"]
        for col in variant_columns:
            if col in results_by_column:
                new_header.append(f"{col}_{results_by_column[col]['Effect']}")
            else:
                new_header.append(col)

        with out_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(new_header)
            writer.writerows(rows)

        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerows(summary_rows)

        self._mark_global_done("annotate_variant_effects")
        return out_path

    def write_effect_summary(self, mut_csv: Path) -> Path:
        out_path = self.output_dir / "effect_category_summary.tsv"
        if self._should_skip_global("effect_category_summary", [out_path]):
            self.logger.info("Skipping effect category summary")
            return out_path

        with mut_csv.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader)

        counts = {"Synonymous": 0, "NonSynonymous": 0, "Unannotated": 0}
        for col in header[1:]:
            if col.endswith("_Synonymous"):
                counts["Synonymous"] += 1
            elif "_NonSynonymous_" in col:
                counts["NonSynonymous"] += 1
            else:
                counts["Unannotated"] += 1

        lines = [
            "Category\tCount",
            f"Synonymous\t{counts['Synonymous']}",
            f"NonSynonymous\t{counts['NonSynonymous']}",
            f"Unannotated\t{counts['Unannotated']}",
        ]
        out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        self._mark_global_done("effect_category_summary")
        return out_path

    def annotate_single_variant(
        self,
        column_name: str,
        cds_map: dict[str, CdsGene],
        cds_sequences: dict[str, str],
    ) -> dict[str, str] | None:
        match = VARIANT_RE.match(column_name)
        if not match:
            return None

        chrom = match.group("chrom")
        pos = int(match.group("pos"))
        ref = match.group("ref")
        alt = match.group("alt")

        if len(ref) != 1 or len(alt) != 1:
            return None

        target_gene = None
        for gene_id, gene in cds_map.items():
            if gene.chrom != chrom:
                continue
            for start, end in gene.cds_regions:
                if start <= pos <= end:
                    target_gene = gene_id
                    break
            if target_gene:
                break

        if target_gene is None:
            return None

        gene = cds_map[target_gene]
        cds_pos = None
        total = 0
        local_ref = ref
        local_alt = alt

        if gene.strand == "+":
            for start, end in gene.cds_regions:
                if start <= pos <= end:
                    cds_pos = total + (pos - start + 1)
                    break
                total += end - start + 1
        else:
            local_ref = self.reverse_complement(ref)
            local_alt = self.reverse_complement(alt)
            for start, end in gene.cds_regions:
                if start <= pos <= end:
                    cds_pos = total + (end - pos + 1)
                    break
                total += end - start + 1

        if cds_pos is None:
            return None

        cds_seq = cds_sequences.get(target_gene, "")
        if cds_pos > len(cds_seq) or cds_pos < 1:
            self.logger.warning(
                "CDS position %s exceeds sequence length for gene %s",
                cds_pos,
                target_gene,
            )
            return None

        current_base = cds_seq[cds_pos - 1]
        if current_base != local_ref:
            self.logger.warning(
                "Reference mismatch for %s at CDS position %s: expected %s, found %s",
                target_gene,
                cds_pos,
                local_ref,
                current_base,
            )
            return None

        codon_start = ((cds_pos - 1) // 3) * 3
        original_codon = cds_seq[codon_start:codon_start + 3]
        if len(original_codon) < 3:
            self.logger.warning("Incomplete codon for %s at CDS position %s", target_gene, cds_pos)
            return None

        pos_in_codon = (cds_pos - 1) % 3
        mutated_codon_list = list(original_codon)
        mutated_codon_list[pos_in_codon] = local_alt
        mutated_codon = "".join(mutated_codon_list)

        original_aa = CODON_TABLE.get(original_codon, "?")
        mutated_aa = CODON_TABLE.get(mutated_codon, "?")
        if original_aa == mutated_aa:
            effect = "Synonymous"
        else:
            effect = f"NonSynonymous_{original_aa}->{mutated_aa}"

        return {
            "Original_Column": column_name,
            "Gene_ID": target_gene,
            "CDS_Position": str(cds_pos),
            "Original_Codon": original_codon,
            "Mutated_Codon": mutated_codon,
            "Effect": effect,
        }

    def parse_gff(self, gff_file: Path) -> dict[str, CdsGene]:
        cds_map: dict[str, CdsGene] = {}
        with gff_file.open("r", encoding="utf-8") as handle:
            for raw in handle:
                if raw.startswith("#"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "CDS":
                    continue
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                parent = self._extract_parent(fields[8])
                if not parent:
                    continue
                parts = parent.split(".")
                gene_id = ".".join(parts[:2]) if len(parts) >= 2 else parent

                if gene_id not in cds_map:
                    cds_map[gene_id] = CdsGene(
                        gene_id=gene_id,
                        chrom=chrom,
                        strand=strand,
                        cds_regions=[],
                    )
                cds_map[gene_id].cds_regions.append((start, end))

        for gene in cds_map.values():
            if gene.strand == "+":
                gene.cds_regions = sorted(gene.cds_regions, key=lambda item: item[0])
            else:
                gene.cds_regions = sorted(gene.cds_regions, key=lambda item: item[0], reverse=True)
        return cds_map

    def parse_fasta(self, fasta_file: Path) -> dict[str, str]:
        sequences: dict[str, list[str]] = {}
        current_id = None
        with fasta_file.open("r", encoding="utf-8") as handle:
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    current_id = line.split()[0][1:]
                    sequences[current_id] = []
                else:
                    if current_id is not None:
                        sequences[current_id].append(line.upper())
        return {key: "".join(value) for key, value in sequences.items()}

    def read_csv_matrix(self, csv_path: Path) -> tuple[list[str], list[list[str]]]:
        with csv_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader)
            rows = [row for row in reader]
        return header, rows

    def reverse_complement(self, seq: str) -> str:
        return seq.upper().translate(COMPLEMENT)[::-1]

    def _extract_parent(self, attributes: str) -> str:
        for item in attributes.split(";"):
            item = item.strip()
            if item.startswith("Parent="):
                return item.split("=", 1)[1]
        return ""

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
        logger = logging.getLogger(f"gwas_step5_{id(self)}")
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="GWAS step5 pipeline: amino-acid effect screening in a target gene region."
    )
    parser.add_argument("--region-vcf", help="Existing regional VCF file")
    parser.add_argument("--vcf", help="Whole-genome VCF or compressed VCF")
    parser.add_argument("--region", help="Target region such as Chr4:49881219-49887065")
    parser.add_argument("--gff", required=True, help="Genome GFF file")
    parser.add_argument("--cds-fasta", required=True, help="CDS FASTA file")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--force", action="store_true", help="Rerun even if outputs already exist")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    pipeline = Step5Pipeline(args)
    pipeline.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
