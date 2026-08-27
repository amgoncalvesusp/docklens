from __future__ import annotations

import pytest

from docklens.parser_mol2 import parse_mol2


def test_gold_mol2_preserves_names_and_first_fitness_score(fixture_path):
    poses = parse_mol2(fixture_path("gold_named_multipose.mol2"))

    assert [pose.molecule_name for pose in poses] == [
        "EOS100588|Actives_100588|mol2|1|dock48",
        "EOS100588|Actives_100588|mol2|1|dock43",
        "EOS100589|Actives_100588|mol2|2|dock2",
    ]
    assert [pose.ligand_id_hint for pose in poses] == [
        "EOS100588",
        "EOS100588",
        "EOS100589",
    ]
    assert [pose.source_pose_label for pose in poses] == [
        "dock48",
        "dock43",
        "dock2",
    ]
    assert [pose.score for pose in poses] == [
        pytest.approx(86.2405),
        pytest.approx(84.125),
        pytest.approx(77.5),
    ]
    assert [pose.score_type for pose in poses] == [
        "Gold.Goldscore.Fitness",
        "Gold.ChemPLP.Fitness",
        "Gold.ASP.Fitness",
    ]


def test_gold_parser_uses_first_valid_fitness_without_treating_gold_score_as_float(
    tmp_path,
):
    path = tmp_path / "scores.mol2"
    path.write_text(
        """@<TRIPOS>MOLECULE
EOS|dock1
1 0 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
1 C1 0 0 0 C.3 1 LIG1 0.0
@<TRIPOS>PROPERTY
> <Gold.Score>
fitness table
> <Gold.Chemscore.Fitness>
11.5
> <Gold.Goldscore.Fitness>
99.0
""",
        encoding="utf-8",
    )

    pose = parse_mol2(path)[0]

    assert pose.score == pytest.approx(11.5)
    assert pose.score_type == "Gold.Chemscore.Fitness"
