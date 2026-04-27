from __future__ import annotations

import json
from pathlib import Path

from padic_grn_viz.engine import load_scores_json, plot_ordering_positions, plot_score_bars


def test_load_scores_json(tmp_path: Path) -> None:
    src = tmp_path / "scores.json"
    src.write_text(json.dumps({"scores": {"mu_E": 1, "mu_A": 2, "mu_I": 3}}), encoding="utf-8")
    payload = load_scores_json(str(src))
    assert payload["scores"]["mu_E"] == 1


def test_plot_outputs(tmp_path: Path) -> None:
    payload = {
        "scores": {"mu_E": 1, "mu_A": 2, "mu_I": 3},
        "ordering": [0, 2, 1],
    }
    dep = plot_score_bars(payload, str(tmp_path / "dep.png"))
    hie = plot_ordering_positions(payload, str(tmp_path / "hie.png"))
    assert dep.exists()
    assert hie.exists()
