"""Exact fast evaluator (p = 2) and genetic algorithm on the three reference networks.

Certified values (exact dynamic program over subsets of genes, mu_dp_exact.js):
  A. thaliana  N = 13: min mu_E = 26,776, attained by exactly four orderings;
                       pi* = [12,3,4,8,5,1,7,2,10,11,9,6,0] has
                       (mu_E, mu_A, mu_I) = (26,776, 53,880, 25,840); max mu_A = 54,560.
  C. elegans   N = 8:  min mu_E = 956, a single ordering [1,0,2,3,6,5,4,7].
  Drosophila   N = 6 (table of this repository): min mu_E = 64, four orderings.
For N = 6 and N = 8 the minimum is recomputed here by exhaustive enumeration.
"""

from __future__ import annotations

import itertools
import json
import random
from pathlib import Path

import pytest

from padic_grn_search import engine
from padic_grn_search.engine import _ReferenceEvaluator, load_reference_core_module, random_ga_search

EX = Path(__file__).resolve().parents[2] / "padic-grn-core" / "examples"
ATH = EX / "transitions_athaliana_N13_p2.csv"
CEL = EX / "transitions_celegans_N8_p2.csv"
DRO = EX / "transitions_drosophila_N6_p2.csv"
PI_STAR = [12, 3, 4, 8, 5, 1, 7, 2, 10, 11, 9, 6, 0]
ATH_OPTIMA = [
    [12, 4, 3, 8, 5, 1, 7, 2, 10, 11, 9, 6, 0],
    [12, 3, 4, 8, 5, 1, 7, 2, 10, 11, 9, 6, 0],
    [12, 4, 3, 8, 1, 5, 7, 2, 10, 11, 9, 6, 0],
    [12, 3, 4, 8, 1, 5, 7, 2, 10, 11, 9, 6, 0],
]

core = load_reference_core_module()


def _calc(csv_path: Path, N: int):
    table = core.read_transition_table_csv(str(csv_path), 2)
    return core.MuCalculator(N=N, p=2, transition_table=table), table


def test_fast_path_is_active_for_complete_boolean_tables() -> None:
    for path, N in ((ATH, 13), (CEL, 8), (DRO, 6)):
        calc, _ = _calc(path, N)
        assert calc._fast is not None


def test_pi_star_and_certified_optima_athaliana() -> None:
    calc, _ = _calc(ATH, 13)
    s = calc.compute_scores_cached(PI_STAR, return_profiles=False)
    assert (s["mu_E"], s["mu_A"], s["mu_I"]) == (26776, 53880, 25840)
    for o in ATH_OPTIMA:
        assert calc.compute_scores_cached(o, return_profiles=False)["mu_E"] == 26776


@pytest.mark.parametrize("path,N,n_ref", [(ATH, 13, 3), (CEL, 8, 150), (DRO, 6, 720)])
def test_fast_equals_reference_and_identity(path: Path, N: int, n_ref: int) -> None:
    calc, table = _calc(path, N)
    rng = random.Random(20261023)
    if N == 6:
        orders = [list(p) for p in itertools.permutations(range(N))]
    else:
        orders = [rng.sample(range(N), N) for _ in range(n_ref)]
    for o in orders[:n_ref]:
        ref = core.compute_scores(2, N, o, table, return_profiles=True)
        fast = calc.compute_scores_cached(o, return_profiles=True)
        for k in ("mu_E", "mu_A", "mu_I", "E", "A", "I"):
            assert fast[k] == ref[k], (o, k)
    for _ in range(500):
        o = rng.sample(range(N), N)
        s = calc.compute_scores_cached(o, return_profiles=False)
        assert s["mu_E"] + s["mu_A"] + s["mu_I"] == N * 2**N


@pytest.mark.parametrize("path,N,best,count", [(CEL, 8, 956, 1), (DRO, 6, 64, 4)])
def test_exhaustive_minimum_small_networks(path: Path, N: int, best: int, count: int) -> None:
    calc, _ = _calc(path, N)
    values = [calc.compute_scores_cached(list(p), return_profiles=False)["mu_E"]
              for p in itertools.permutations(range(N))]
    assert min(values) == best
    assert values.count(best) == count


@pytest.mark.parametrize("path,N", [(ATH, 13), (CEL, 8)])
def test_incremental_moves_equal_full_evaluation(path: Path, N: int) -> None:
    calc, table = _calc(path, N)
    full, _ = _calc(path, N)
    rng = random.Random(7)
    for k in range(300):
        A = rng.sample(range(N), N)
        i, j = sorted(rng.sample(range(N), 2))
        B = list(A)
        B[i], B[j] = B[j], B[i]
        assert calc.compute_scores_swap(A, i, j, return_profiles=True) == full.compute_scores_cached(B, return_profiles=True)
        C = list(A)
        C.insert(j, C.pop(i))
        assert calc.compute_scores_neighbor(A, C, return_profiles=True) == full.compute_scores_cached(C, return_profiles=True)


def test_ga_finds_certified_minimum_athaliana_default_mode() -> None:
    out = random_ga_search(str(ATH), 13, 2, generations=50, population=30, seed=1)
    assert out["best_fitness"] == 26776
    assert out["best_ordering"] in ATH_OPTIMA
    s = out["best_scores"]
    assert s["mu_E"] + s["mu_A"] + s["mu_I"] == 13 * 2**13


def test_ga_memetic_mode_finds_minimum_small_budget() -> None:
    for seed in (1, 2, 3):
        out = random_ga_search(str(ATH), 13, 2, generations=10, population=20, seed=seed, local_search=True)
        assert out["best_fitness"] == 26776
        assert out["method"] == "ga_pmx_tournament_elite_memetic"


def test_ga_maximize_mu_a_memetic() -> None:
    out = random_ga_search(str(ATH), 13, 2, objective="maximize_mu_A", generations=20, population=20,
                           seed=1, local_search=True)
    assert out["best_scores"]["mu_A"] == 54560


def test_ga_reproducible_by_seed() -> None:
    a = random_ga_search(str(CEL), 8, 2, generations=30, population=20, seed=42)
    b = random_ga_search(str(CEL), 8, 2, generations=30, population=20, seed=42)
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)


def test_ga_same_trajectory_with_and_without_fast_path(monkeypatch: pytest.MonkeyPatch) -> None:
    """The fast path changes no random draw and no score, so the whole run is identical."""
    fast = random_ga_search(str(CEL), 8, 2, generations=25, population=16, seed=5)
    monkeypatch.setattr(core, "table_supports_fast_path", lambda *a, **k: False)
    monkeypatch.setattr(engine, "load_reference_core_module", lambda: core)
    ev = _ReferenceEvaluator(str(CEL), 8, 2)
    assert ev.calculator._fast is None
    slow = random_ga_search(str(CEL), 8, 2, generations=25, population=16, seed=5)
    assert json.dumps(fast, sort_keys=True) == json.dumps(slow, sort_keys=True)
