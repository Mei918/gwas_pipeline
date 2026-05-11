from __future__ import annotations

import argparse
import json
import logging
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional


FASTQ_PATTERN = re.compile(r"^(?P<sample>.+)\.R(?P<read>[12])\.fastq\.gz$")


@dataclass
class SamplePair:
    sample_id: str
    read1: Path
    read2: Path


@dataclass
class PipelineConfig:
    fastq_dir: Path
    reference_fasta: Path
    output_dir: Path
    fastp_threads: int = 8
    bwa_threads: int = 32
    samtools_sort_threads: int = 8
    maf: float = 0.05
    max_missing: float = 0.8
    force: bool = False


class Step1Pipeline:
    def __init__(self, config: PipelineConfig) -> None:
        self.config = config
        self.dirs = self._ensure_dirs()
        self.logger = self._build_logger()
        self.state_path = self.dirs["state"] / "pipeline_state.json"
        self.state = self._load_state()

    def run(self) -> None:
        self.logger.info("Starting GWAS step1 pipeline")
        samples = self.scan_samples()
        if not samples:
            raise RuntimeError(f"No paired FASTQ files found in {self.config.fastq_dir}")

        self.logger.info("Discovered %d samples", len(samples))
        self.build_reference_index()

        for sample in samples:
            self.run_fastp(sample)
            self.align_sample(sample)
            self.call_sample_gvcf(sample)

        gvcf_files = [self.sample_gvcf_path(sample.sample_id) for sample in samples]
        self.combine_gvcfs(gvcf_files)
        self.genotype_cohort()
        self.filter_vcf()
        self.make_plink()
        self.logger.info("Pipeline finished successfully")

    def scan_samples(self) -> list[SamplePair]:
        found: dict[str, dict[str, Path]] = {}
        for fastq in sorted(self.config.fastq_dir.glob("*.fastq.gz")):
            match = FASTQ_PATTERN.match(fastq.name)
            if not match:
                self.logger.warning("Skipping unmatched file: %s", fastq.name)
                continue
            sample_id = match.group("sample")
            read_key = match.group("read")
            if sample_id not in found:
                found[sample_id] = {}
            found[sample_id][read_key] = fastq

        pairs = []
        for sample_id, reads in sorted(found.items()):
            read1 = reads.get("1")
            read2 = reads.get("2")
            if read1 and read2:
                pairs.append(SamplePair(sample_id=sample_id, read1=read1, read2=read2))
            else:
                self.logger.warning("Skipping incomplete pair for sample %s", sample_id)
        return pairs

    def build_reference_index(self) -> None:
        ref = self.config.reference_fasta
        bwa_index_files = [Path(f"{ref}.{suffix}") for suffix in ["amb", "ann", "bwt", "pac", "sa"]]
        fai_file = Path(f"{ref}.fai")
        dict_file = ref.with_suffix(".dict")

        if self._should_skip_global("reference_index", bwa_index_files + [fai_file, dict_file]):
            self.logger.info("Skipping reference indexing")
            return

        self._cleanup_paths_if_force(bwa_index_files + [fai_file, dict_file])
        self._run_cmd(f"bwa index {self._q(ref)}")
        self._run_cmd(f"samtools faidx {self._q(ref)}")
        self._run_cmd(
            f"gatk CreateSequenceDictionary -R {self._q(ref)} -O {self._q(dict_file)}"
        )
        self._mark_global_done("reference_index")

    def run_fastp(self, sample: SamplePair) -> None:
        clean_r1 = self.clean_read1_path(sample.sample_id)
        clean_r2 = self.clean_read2_path(sample.sample_id)
        html = self.dirs["qc"] / f"{sample.sample_id}.html"
        json_report = self.dirs["qc"] / f"{sample.sample_id}.json"

        if self._should_skip_sample(sample.sample_id, "fastp", [clean_r1, clean_r2, html, json_report]):
            self.logger.info("Skipping fastp for %s", sample.sample_id)
            return

        cmd = (
            f"fastp "
            f"-i {self._q(sample.read1)} "
            f"-I {self._q(sample.read2)} "
            f"-o {self._q(clean_r1)} "
            f"-O {self._q(clean_r2)} "
            f"-q 20 -u 30 -n 5 -l 50 "
            f"-w {self.config.fastp_threads} "
            f"-h {self._q(html)} "
            f"-j {self._q(json_report)}"
        )
        self._run_cmd(cmd)
        self._mark_sample_done(sample.sample_id, "fastp")

    def align_sample(self, sample: SamplePair) -> None:
        bam_file = self.sample_bam_path(sample.sample_id)
        bai_file = Path(f"{bam_file}.bai")
        clean_r1 = self.clean_read1_path(sample.sample_id)
        clean_r2 = self.clean_read2_path(sample.sample_id)

        if self._should_skip_sample(sample.sample_id, "align", [bam_file, bai_file]):
            self.logger.info("Skipping alignment for %s", sample.sample_id)
            return

        cmd = (
            f"bwa mem -t {self.config.bwa_threads} "
            f"{self._q(self.config.reference_fasta)} "
            f"{self._q(clean_r1)} "
            f"{self._q(clean_r2)} | "
            f"samtools sort -@ {self.config.samtools_sort_threads} -o {self._q(bam_file)}"
        )
        self._run_cmd(cmd)
        rg_bam = bam_file.with_suffix(".rg.bam")
        self._run_cmd(
            f"gatk AddOrReplaceReadGroups "
            f"-I {self._q(bam_file)} "
            f"-O {self._q(rg_bam)} "
            f"-RGID {sample.sample_id} "
            f"-RGLB lib1 "
            f"-RGPL ILLUMINA "
            f"-RGPU unit1 "
            f"-RGSM {sample.sample_id}"
        )
        rg_bam.replace(bam_file)
        self._run_cmd(f"samtools index {self._q(bam_file)}")
        self._mark_sample_done(sample.sample_id, "align")

    def call_sample_gvcf(self, sample: SamplePair) -> None:
        bam_file = self.sample_bam_path(sample.sample_id)
        gvcf_file = self.sample_gvcf_path(sample.sample_id)

        if self._should_skip_sample(sample.sample_id, "gvcf", [gvcf_file]):
            self.logger.info("Skipping GVCF calling for %s", sample.sample_id)
            return

        cmd = (
            f"gatk HaplotypeCaller "
            f"-R {self._q(self.config.reference_fasta)} "
            f"-I {self._q(bam_file)} "
            f"-O {self._q(gvcf_file)} "
            f"-ERC GVCF"
        )
        self._run_cmd(cmd)
        self._mark_sample_done(sample.sample_id, "gvcf")

    def combine_gvcfs(self, gvcf_files: list[Path]) -> None:
        cohort_gvcf = self.cohort_gvcf_path()
        if self._should_skip_global("combine_gvcf", [cohort_gvcf]):
            self.logger.info("Skipping GVCF combine")
            return

        missing = [path for path in gvcf_files if not path.exists()]
        if missing:
            raise RuntimeError(f"Missing per-sample GVCF files: {missing[:5]}")

        variants = " ".join(f"-V {self._q(path)}" for path in gvcf_files)
        cmd = (
            f"gatk CombineGVCFs "
            f"-R {self._q(self.config.reference_fasta)} "
            f"{variants} "
            f"-O {self._q(cohort_gvcf)}"
        )
        self._run_cmd(cmd)
        self._mark_global_done("combine_gvcf")

    def genotype_cohort(self) -> None:
        raw_vcf = self.raw_vcf_path()
        if self._should_skip_global("genotype_gvcf", [raw_vcf]):
            self.logger.info("Skipping cohort genotyping")
            return

        cmd = (
            f"gatk GenotypeGVCFs "
            f"-R {self._q(self.config.reference_fasta)} "
            f"-V {self._q(self.cohort_gvcf_path())} "
            f"-O {self._q(raw_vcf)}"
        )
        self._run_cmd(cmd)
        self._mark_global_done("genotype_gvcf")

    def filter_vcf(self) -> None:
        filtered_vcf = self.filtered_vcf_path()
        if self._should_skip_global("filter_vcf", [filtered_vcf]):
            self.logger.info("Skipping VCF filtering")
            return

        out_prefix = self.dirs["vcf"] / "cohort.gwas"
        cmd = (
            f"vcftools --gzvcf {self._q(self.raw_vcf_path())} "
            f"--maf {self.config.maf} "
            f"--max-missing {self.config.max_missing} "
            f"--recode --out {self._q(out_prefix)}"
        )
        self._run_cmd(cmd)
        self._mark_global_done("filter_vcf")

    def make_plink(self) -> None:
        bed = self.dirs["plink"] / "cohort.bed"
        bim = self.dirs["plink"] / "cohort.bim"
        fam = self.dirs["plink"] / "cohort.fam"
        if self._should_skip_global("make_plink", [bed, bim, fam]):
            self.logger.info("Skipping PLINK conversion")
            return

        if self._count_vcf_variants(self.filtered_vcf_path()) == 0:
            self.logger.warning("Filtered VCF has no variants; skipping PLINK conversion")
            self._mark_global_done("make_plink")
            return

        out_prefix = self.dirs["plink"] / "cohort"
        cmd = (
            f"plink --vcf {self._q(self.filtered_vcf_path())} "
            f"--make-bed --out {self._q(out_prefix)}"
        )
        self._run_cmd(cmd)
        self._mark_global_done("make_plink")

    def clean_read1_path(self, sample_id: str) -> Path:
        return self.dirs["clean"] / f"{sample_id}.R1.clean.fastq.gz"

    def clean_read2_path(self, sample_id: str) -> Path:
        return self.dirs["clean"] / f"{sample_id}.R2.clean.fastq.gz"

    def sample_bam_path(self, sample_id: str) -> Path:
        return self.dirs["bam"] / f"{sample_id}.bam"

    def sample_gvcf_path(self, sample_id: str) -> Path:
        return self.dirs["gvcf"] / f"{sample_id}.g.vcf.gz"

    def cohort_gvcf_path(self) -> Path:
        return self.dirs["vcf"] / "cohort.g.vcf.gz"

    def raw_vcf_path(self) -> Path:
        return self.dirs["vcf"] / "cohort.raw.vcf.gz"

    def filtered_vcf_path(self) -> Path:
        return self.dirs["vcf"] / "cohort.gwas.recode.vcf"

    def _ensure_dirs(self) -> dict[str, Path]:
        base = self.config.output_dir
        dirs = {
            "base": base,
            "logs": base / "logs",
            "state": base / "state",
            "clean": base / "clean",
            "qc": base / "qc",
            "bam": base / "bam",
            "gvcf": base / "gvcf",
            "vcf": base / "vcf",
            "plink": base / "plink",
        }
        for path in dirs.values():
            path.mkdir(parents=True, exist_ok=True)
        return dirs

    def _build_logger(self) -> logging.Logger:
        logger = logging.getLogger(f"gwas_step1_{id(self)}")
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
        return {"global": {}, "samples": {}}

    def _save_state(self) -> None:
        self.state_path.write_text(json.dumps(self.state, indent=2, sort_keys=True))

    def _mark_sample_done(self, sample_id: str, step: str) -> None:
        self.state["samples"].setdefault(sample_id, {})
        self.state["samples"][sample_id][step] = "done"
        self._save_state()

    def _mark_global_done(self, step: str) -> None:
        self.state["global"][step] = "done"
        self._save_state()

    def _should_skip_sample(self, sample_id: str, step: str, outputs: list[Path]) -> bool:
        if self.config.force:
            return False
        state_done = self.state["samples"].get(sample_id, {}).get(step) == "done"
        files_done = all(path.exists() for path in outputs)
        return state_done or files_done

    def _should_skip_global(self, step: str, outputs: list[Path]) -> bool:
        if self.config.force:
            return False
        state_done = self.state["global"].get(step) == "done"
        files_done = all(path.exists() for path in outputs)
        return state_done or files_done

    def _cleanup_paths_if_force(self, paths: list[Path]) -> None:
        if not self.config.force:
            return
        for path in paths:
            if path.exists():
                path.unlink()

    def _count_vcf_variants(self, path: Path) -> int:
        if not path.exists():
            return 0
        count = 0
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith("#"):
                    count += 1
        return count

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
        description="GWAS step1 pipeline: fastq.gz to filtered VCF and PLINK."
    )
    parser.add_argument("--fastq-dir", required=True, help="Directory containing paired FASTQ files")
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--fastp-threads", type=int, default=8, help="Threads for fastp")
    parser.add_argument("--bwa-threads", type=int, default=32, help="Threads for bwa mem")
    parser.add_argument("--sort-threads", type=int, default=8, help="Threads for samtools sort")
    parser.add_argument("--maf", type=float, default=0.05, help="MAF threshold for vcftools")
    parser.add_argument("--max-missing", type=float, default=0.8, help="Max missing threshold for vcftools")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rerun steps even if outputs or state already exist",
    )
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    config = PipelineConfig(
        fastq_dir=Path(args.fastq_dir).resolve(),
        reference_fasta=Path(args.reference).resolve(),
        output_dir=Path(args.output_dir).resolve(),
        fastp_threads=args.fastp_threads,
        bwa_threads=args.bwa_threads,
        samtools_sort_threads=args.sort_threads,
        maf=args.maf,
        max_missing=args.max_missing,
        force=args.force,
    )
    pipeline = Step1Pipeline(config)
    pipeline.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
