from __future__ import annotations

from pathlib import Path

from padic_grn_search.engine import branch_and_bound_search, exhaustive_search, random_ga_search


def _write_identity_csv(path: Path, p: int, n: int) -> None:
    total = p**n
    rows = ["config_from,config_to"]
    rows.extend(f"{i},{i}" for i in range(total))
    path.write_text("\n".join(rows), encoding="utf-8")


def test_random_ga_search_returns_result(tmp_path: Path) -> None:
    csv_path = tmp_path / "transitions.csv"
    _write_identity_csv(csv_path, p=2, n=2)
    out = random_ga_search(
        transitions_csv=str(csv_path),
        n=2,
        p=2,
        generations=2,
        population=4,
        seed=7,
    )
    assert out["best_ordering"] in ([0, 1], [1, 0])
    assert set(out["best_scores"].keys()) == {"mu_E", "mu_A", "mu_I"}


def test_exhaustive_search_truncated(tmp_path: Path) -> None:
    csv_path = tmp_path / "transitions.csv"
    _write_identity_csv(csv_path, p=2, n=3)
    out = exhaustive_search(
        transitions_csv=str(csv_path),
        n=3,
        p=2,
        max_evaluations=2,
    )
    assert out["evaluations"] == 2
    assert out["truncated"] is True


def test_branch_and_bound_search_returns_result(tmp_path: Path) -> None:
    csv_path = tmp_path / "transitions.csv"
    _write_identity_csv(csv_path, p=2, n=2)
    out = branch_and_bound_search(
        transitions_csv=str(csv_path),
        n=2,
        p=2,
        objective="minimize_mu_E",
    )
    assert out["best_ordering"] in ([0, 1], [1, 0])
    assert out["best_scores"]["mu_E"] == 0
