"""Step registry for the pure GWAS pipeline package."""

from __future__ import annotations

from importlib import import_module
from typing import Callable, Optional


StepMain = Callable[[Optional[list[str]]], int]

STEP_MODULES: dict[str, str] = {
    "step1": "gwas_pipeline.steps.step1_fastq_to_vcf",
    "step2": "gwas_pipeline.steps.step2_gwas_analysis",
    "step3": "gwas_pipeline.steps.step3_candidate_gene_extraction",
    "step4": "gwas_pipeline.steps.step4_local_manhattan_ld",
    "step5": "gwas_pipeline.steps.step5_amino_acid_mutation",
    "step6": "gwas_pipeline.steps.step6_haplotype_analysis",
    "step7": "gwas_pipeline.steps.step7_haplotype_network",
    "step8": "gwas_pipeline.steps.step8_pi_tajima",
    "step9": "gwas_pipeline.steps.step9_geo_visualization",
}


def get_step_main(step_name: str) -> StepMain:
    try:
        module_name = STEP_MODULES[step_name]
    except KeyError as exc:
        raise KeyError(f"Unknown GWAS step: {step_name}") from exc

    module = import_module(module_name)
    return module.main


__all__ = ["STEP_MODULES", "StepMain", "get_step_main"]
