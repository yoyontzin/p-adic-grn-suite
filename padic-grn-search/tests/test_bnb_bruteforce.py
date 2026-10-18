"""Branch and bound contra fuerza bruta en redes aleatorias chicas (niveles n = 0..N-1)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from padic_grn_search.engine import branch_and_bound_search, evaluate_ordering, exhaustive_search

_BF_PATH = Path(__file__).resolve().parents[2] / "padic-grn-core" / "tests" / "_bruteforce.py"
_spec = importlib.util.spec_from_file_location("_bruteforce_core", _BF_PATH)
bf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bf)


def _write_csv(path: Path, table: dict[int, int]) -> None:
    rows = ["config_from,config_to"] + [f"{x},{y}" for x, y in sorted(table.items())]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


@pytest.mark.parametrize("N,p,seed,n_const", bf.CASES)
def test_bnb_optima_match_bruteforce(tmp_path: Path, N: int, p: int, seed: int, n_const: int) -> None:
    table = bf.random_network(N, p, seed, n_const)
    csv_path = tmp_path / "t.csv"
    _write_csv(csv_path, table)
    best_E, best_A = bf.optima(table, N, p)

    out_min = branch_and_bound_search(str(csv_path), N, p, objective="minimize_mu_E")
    assert out_min["method"] == "branch_and_bound_partial_muE"
    assert out_min["best_fitness"] == best_E
    assert bf.scores(table, N, p, out_min["best_ordering"])[0] == best_E

    out_max = branch_and_bound_search(str(csv_path), N, p, objective="maximize_mu_A")
    assert out_max["method"] == "branch_and_bound_upper_muA"
    assert out_max["best_fitness"] == -best_A
    assert bf.scores(table, N, p, out_max["best_ordering"])[1] == best_A
    assert out_max["best_scores"]["mu_A"] == best_A

    for out in (out_min, out_max):
        s = out["best_scores"]
        assert s["mu_E"] + s["mu_A"] + s["mu_I"] == N * p**N


@pytest.mark.parametrize("N,p,seed,n_const", [c for c in bf.CASES if c[0] <= 5])
def test_exhaustive_matches_bruteforce(tmp_path: Path, N: int, p: int, seed: int, n_const: int) -> None:
    table = bf.random_network(N, p, seed, n_const)
    csv_path = tmp_path / "t.csv"
    _write_csv(csv_path, table)
    best_E, best_A = bf.optima(table, N, p)
    assert exhaustive_search(str(csv_path), N, p, objective="minimize_mu_E")["best_fitness"] == best_E
    assert exhaustive_search(str(csv_path), N, p, objective="maximize_mu_A")["best_fitness"] == -best_A


def test_pi_star_athaliana() -> None:
    csv_path = Path(__file__).resolve().parents[2] / "padic-grn-core" / "examples" / "transitions_athaliana_N13_p2.csv"
    scores = evaluate_ordering(str(csv_path), 13, 2, [12, 3, 4, 8, 5, 1, 7, 2, 10, 11, 9, 6, 0])
    assert scores == {"mu_E": 26776, "mu_A": 53880, "mu_I": 25840}
