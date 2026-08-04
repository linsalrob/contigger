"""Optional minimap2 synthetic integration tests."""

from pathlib import Path

import pytest

from contigger.aligners.minimap2 import Minimap2Aligner
from contigger.config import build_run_config
from contigger.relationships import classify_pair, group_ordered_pairs
from contigger.synthetic import synthetic_cases
from contigger.utilities.subprocesses import find_executable

MINIMAP2 = find_executable("minimap2")


@pytest.mark.integration
@pytest.mark.skipif(MINIMAP2 is None, reason="minimap2 is not installed")
@pytest.mark.parametrize("preset", ["asm5", "asm20"])
def test_minimap2_presets_on_synthetic_truth(tmp_path: Path, preset: str) -> None:
    config = build_run_config(min_overlap=100, min_containment=100, end_tolerance=10)
    simple_cases = synthetic_cases()[:8]
    for case in simple_cases:
        target = tmp_path / f"{case.name}.target.fa"
        query = tmp_path / f"{case.name}.query.fa"
        target.write_text(f">target\n{case.target}\n", encoding="ascii")
        query.write_text(f">query\n{case.query}\n", encoding="ascii")
        aligner = Minimap2Aligner(executable=MINIMAP2, preset=preset)
        groups = list(group_ordered_pairs(aligner.align_paths(target, query)))
        assert groups, f"{preset} missed fixture {case.name}"
        decision = classify_pair(groups[0], config)
        assert decision.relationship.relationship_type is case.truth
        assert aligner.tool_version
        assert aligner.last_command is not None


def test_minimap2_default_preset_remains_asm20() -> None:
    aligner = Minimap2Aligner(executable=Path("/usr/bin/minimap2"))
    assert aligner.preset == "asm20"
    assert "asm20" in aligner.command_for_paths(Path("target.fa"), Path("query.fa"))
