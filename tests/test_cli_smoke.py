from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_cmd(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args,
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


def test_module_help() -> None:
    result = run_cmd(sys.executable, "-m", "gwas_pipeline", "--help")
    assert result.returncode == 0
    assert "gwas-pipeline" in result.stdout
    assert "step9" in result.stdout


def test_doctor_help() -> None:
    result = run_cmd(sys.executable, "-m", "gwas_pipeline", "doctor", "--help")
    assert result.returncode == 0
    assert "plink_env" in result.stdout


def test_step_help() -> None:
    result = run_cmd(sys.executable, "-m", "gwas_pipeline", "step8", "--help")
    assert result.returncode == 0
    assert "--vcf" in result.stdout


def test_legacy_wrapper_help() -> None:
    wrapper = ROOT / "gwas_step1_project" / "step1_fastq_to_vcf.py"
    result = run_cmd(sys.executable, str(wrapper), "--help")
    assert result.returncode == 0
    assert "--fastq-dir" in result.stdout
