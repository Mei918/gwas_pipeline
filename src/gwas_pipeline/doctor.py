from __future__ import annotations

import argparse
import importlib
import shutil
import sys
from dataclasses import dataclass
from typing import Optional


@dataclass
class CheckResult:
    name: str
    status: str
    detail: str


TOOL_GROUPS = {
    "plink_env": ["plink", "bcftools", "vcftools", "bwa", "samtools", "fastp", "gatk", "gcta64"],
    "omicverse": [],
}

PYTHON_GROUPS = {
    "omicverse": ["pandas", "numpy", "matplotlib", "networkx", "yaml"],
    "base": [],
}


def run_doctor(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="gwas-pipeline doctor",
        description="Check whether server-side GWAS package dependencies are available.",
    )
    parser.add_argument(
        "--profile",
        choices=["all", "plink_env", "omicverse"],
        default="all",
        help="Subset of checks to run.",
    )
    args = parser.parse_args(argv)

    profiles = [args.profile] if args.profile != "all" else ["plink_env", "omicverse"]
    results: list[CheckResult] = []

    for profile in profiles:
        for tool in TOOL_GROUPS.get(profile, []):
            path = shutil.which(tool)
            results.append(
                CheckResult(
                    name=f"{profile}:tool:{tool}",
                    status="OK" if path else "MISSING",
                    detail=path or "not found on PATH",
                )
            )

        for module_name in PYTHON_GROUPS.get(profile, []):
            try:
                mod = importlib.import_module(module_name)
                detail = getattr(mod, "__file__", "imported")
                results.append(CheckResult(name=f"{profile}:python:{module_name}", status="OK", detail=str(detail)))
            except Exception as exc:
                results.append(CheckResult(name=f"{profile}:python:{module_name}", status="MISSING", detail=str(exc)))

    ok = 0
    missing = 0
    for result in results:
        if result.status == "OK":
            ok += 1
        else:
            missing += 1
        print(f"{result.status}\t{result.name}\t{result.detail}")

    print(f"summary\tok={ok}\tmissing={missing}")
    return 0 if missing == 0 else 1


if __name__ == "__main__":
    raise SystemExit(run_doctor())
