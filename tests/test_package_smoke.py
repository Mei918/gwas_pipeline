from __future__ import annotations

from gwas_pipeline import __version__
from gwas_pipeline.steps import STEP_MODULES, get_step_main


def test_version_is_present() -> None:
    assert isinstance(__version__, str)
    assert __version__


def test_step_registry_contains_all_steps() -> None:
    assert list(STEP_MODULES) == [
        "step1",
        "step2",
        "step3",
        "step4",
        "step5",
        "step6",
        "step7",
        "step8",
        "step9",
    ]


def test_each_step_resolves_to_callable_main() -> None:
    for step_name in STEP_MODULES:
        step_main = get_step_main(step_name)
        assert callable(step_main), step_name
