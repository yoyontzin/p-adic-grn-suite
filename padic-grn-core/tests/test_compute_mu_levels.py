"""Motor de padic-grn-core contra fuerza bruta: niveles n = 0..N-1 y cotas parciales."""

from __future__ import annotations

import importlib.util
import itertools
import random
import shutil
import subprocess
from pathlib import Path

import pytest

import _bruteforce as bf
from padic_grn_core import compute_mu

ROOT = Path(__file__).resolve().parents[1]


def _orderings(N: int, rng: random.Random, limit: int = 60) -> list[list[int]]:
    perms = [list(A) for A in itertools.permutations(range(N))]
    return perms if len(perms) <= limit else rng.sample(perms, limit)


@pytest.mark.parametrize("N,p,seed,n_const", bf.CASES)
def test_scores_match_bruteforce(N: int, p: int, seed: int, n_const: int) -> None:
    table = bf.random_network(N, p, seed, n_const)
    assert len(bf.constant_genes(table, N, p)) >= n_const
    rng = random.Random(seed)
    for A in _orderings(N, rng):
        out = compute_mu.compute_scores(p, N, A, table, return_profiles=True)
        cached = compute_mu.MuCalculator(N, p, table).compute_scores_cached(A)
        counts = bf.level_counts(table, N, p, A)
        assert (out["mu_E"], out["mu_A"], out["mu_I"]) == bf.scores(table, N, p, A)
        assert {k: cached[k] for k in ("mu_E", "mu_A", "mu_I")} == {k: out[k] for k in ("mu_E", "mu_A", "mu_I")}
        assert sorted(out["E"]) == list(range(N))
        for n in range(N):
            assert (out["E"][n], out["A"][n], out["I"][n]) == counts[n]
            assert sum(counts[n]) == p**n
        assert out["mu_E"] + out["mu_A"] + out["mu_I"] == N * p**N


@pytest.mark.parametrize("N,p,seed,n_const", [c for c in bf.CASES if c[3] >= 1])
def test_root_contracts_iff_first_gene_constant(N: int, p: int, seed: int, n_const: int) -> None:
    table = bf.random_network(N, p, seed, n_const)
    const = bf.constant_genes(table, N, p)
    for A in _orderings(N, random.Random(seed), limit=30):
        out = compute_mu.compute_scores(p, N, A, table)
        assert out["E"][0] == 0
        if A[0] in const:
            assert (out["A"][0], out["I"][0]) == (1, 0)
        else:
            assert (out["A"][0], out["I"][0]) == (0, 1)


@pytest.mark.parametrize("N,p,seed,n_const", bf.CASES)
def test_partial_bounds(N: int, p: int, seed: int, n_const: int) -> None:
    table = bf.random_network(N, p, seed, n_const)
    rng = random.Random(seed)
    for A in _orderings(N, rng, limit=12):
        for k in range(-1, N):
            prefix = A[: k + 1]
            counts = compute_mu.partial_level_counts(prefix, N, p, table)
            assert sorted(counts) == list(range(k + 1))
            lower_E = compute_mu.partial_mu_E(prefix, N, p, table)
            lower_A = compute_mu.partial_mu_A(prefix, N, p, table)
            upper_A = compute_mu.mu_A_upper_bound(prefix, N, p, table)
            assert upper_A == lower_A + (N - 1 - k) * p**N
            rest = [g for g in range(N) if g not in prefix]
            completions = [prefix + list(t) for t in itertools.permutations(rest)]
            for full in rng.sample(completions, min(8, len(completions))):
                full_counts = bf.level_counts(table, N, p, full)
                assert all(counts[n] == full_counts[n] for n in counts)
                mu_E, mu_A, _ = bf.scores(table, N, p, full)
                assert lower_E <= mu_E
                assert lower_A <= mu_A <= upper_A
            if k == N - 1:
                mu_E, mu_A, _ = bf.scores(table, N, p, A)
                assert (lower_E, lower_A, upper_A) == (mu_E, mu_A, mu_A)


def test_pi_star_athaliana() -> None:
    table = compute_mu.read_transition_table_csv(str(ROOT / "examples" / "transitions_athaliana_N13_p2.csv"), 2)
    assert len(table) == 2**13
    out = compute_mu.compute_scores(2, 13, [12, 3, 4, 8, 5, 1, 7, 2, 10, 11, 9, 6, 0], table)
    assert (out["mu_E"], out["mu_A"], out["mu_I"]) == (26776, 53880, 25840)


@pytest.mark.skipif(shutil.which("node") is None, reason="node no disponible")
def test_pi_star_athaliana_dp_js() -> None:
    res = subprocess.run(
        ["node", str(ROOT / "tools" / "mu_dp_exact.js"), str(ROOT / "examples" / "transitions_athaliana_N13_p2.csv"),
         "--eval", "12,3,4,8,5,1,7,2,10,11,9,6,0"],
        capture_output=True, text=True, check=True, timeout=60,
    )
    assert "(mu_E, mu_A, mu_I) = (26776, 53880, 25840)" in res.stdout


def test_find_optimal_orderings_tool() -> None:
    spec = importlib.util.spec_from_file_location("find_optimal_orderings", ROOT / "tools" / "find_optimal_orderings.py")
    tool = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tool)
    for N, p, seed, n_const in bf.CASES[:6]:
        table = bf.random_network(N, p, seed, n_const)
        tt = [table[x] for x in range(p**N)]
        for A in _orderings(N, random.Random(seed), limit=10):
            assert tool.compute_mu(tt, N, p, A) == bf.scores(table, N, p, A)
    assert tool.compute_mu(tool.TT_DROSO, 6, 2, [0, 5, 3, 4, 2, 1]) == (64, 156, 164)
