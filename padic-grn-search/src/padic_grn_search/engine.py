"""Motor de búsqueda para padic-grn-search.

Incluye:
- Algoritmo genético con selección por torneo, PMX, mutación y élite.
- Branch-and-Bound / backtracking exacto para minimizar mu_E o maximizar mu_A.
- Búsqueda exhaustiva (baseline y verificación).
"""

from __future__ import annotations

import importlib.util
import itertools
import random
from pathlib import Path
from typing import Any


def _vendored_core_path() -> Path:
    # engine.py -> padic_grn_search -> src -> padic-grn-search -> raíz de la suite
    return Path(__file__).resolve().parents[3] / "padic-grn-core" / "src" / "padic_grn_core" / "compute_mu.py"


def load_reference_core_module() -> Any:
    """Motor de referencia de padic-grn-core (niveles n = 0..N-1)."""
    try:
        from padic_grn_core import compute_mu

        return compute_mu
    except ImportError:
        pass
    module_path = _vendored_core_path()
    if not module_path.exists():
        raise FileNotFoundError(f"No se encontró el motor de padic-grn-core: {module_path}")
    spec = importlib.util.spec_from_file_location("padic_grn_core_compute_mu", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("No se pudo crear spec del módulo core de referencia.")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _ReferenceEvaluator:
    def __init__(self, transitions_csv: str, n: int, p: int):
        self.core = load_reference_core_module()
        self.transitions_csv = transitions_csv
        self.n = n
        self.p = p
        self.table = self.core.read_transition_table_csv(transitions_csv, p)
        self.calculator = None
        if hasattr(self.core, "MuCalculator"):
            self.calculator = self.core.MuCalculator(N=n, p=p, transition_table=self.table)

    def score(self, ordering: list[int]) -> dict[str, int]:
        if self.calculator is not None:
            scores = self.calculator.compute_scores_cached(ordering, return_profiles=False)
        else:
            scores = self.core.compute_scores(self.p, self.n, ordering, self.table, return_profiles=False)
        return {"mu_E": scores["mu_E"], "mu_A": scores["mu_A"], "mu_I": scores["mu_I"]}

    def partial_mu_e(self, prefix: list[int]) -> int:
        """Cota inferior de mu_E: niveles n = 0..k fijados por el prefijo."""
        return int(self.core.partial_mu_E(prefix, self.n, self.p, self.table))

    def mu_a_upper_bound(self, prefix: list[int]) -> int:
        """Cota superior de mu_A: partial_mu_A + (N-1-k) p^N."""
        return int(self.core.mu_A_upper_bound(prefix, self.n, self.p, self.table))


def evaluate_ordering(
    transitions_csv: str,
    n: int,
    p: int,
    ordering: list[int],
) -> dict[str, int]:
    evaluator = _ReferenceEvaluator(transitions_csv, n, p)
    return evaluator.score(ordering)


def _fitness(scores: dict[str, int], objective: str) -> int:
    return scores["mu_E"] if objective == "minimize_mu_E" else -scores["mu_A"]


def _pmx_crossover(parent1: list[int], parent2: list[int], rng: random.Random) -> tuple[list[int], list[int]]:
    n = len(parent1)
    c1, c2 = sorted(rng.sample(range(n), 2))

    def build_child(p1: list[int], p2: list[int]) -> list[int]:
        child = [-1] * n
        child[c1 : c2 + 1] = p1[c1 : c2 + 1]
        for i in range(c1, c2 + 1):
            gene = p2[i]
            if gene in child:
                continue
            pos = i
            while c1 <= pos <= c2:
                pos = p2.index(p1[pos])
            child[pos] = gene
        for i in range(n):
            if child[i] == -1:
                child[i] = p2[i]
        return child

    return build_child(parent1, parent2), build_child(parent2, parent1)


def _mutate_swap(ordering: list[int], rng: random.Random, mutation_rate: float) -> list[int]:
    out = list(ordering)
    for i in range(len(out)):
        if rng.random() < mutation_rate:
            j = rng.randrange(len(out))
            out[i], out[j] = out[j], out[i]
    return out


def _tournament_select(
    population: list[list[int]],
    fitnesses: list[int],
    rng: random.Random,
    k: int = 3,
) -> list[int]:
    idxs = [rng.randrange(len(population)) for _ in range(k)]
    best_idx = min(idxs, key=lambda idx: fitnesses[idx])
    return list(population[best_idx])


def random_ga_search(
    transitions_csv: str,
    n: int,
    p: int,
    objective: str = "minimize_mu_E",
    generations: int = 50,
    population: int = 30,
    seed: int | None = None,
) -> dict[str, Any]:
    rng = random.Random(seed)
    evaluator = _ReferenceEvaluator(transitions_csv, n, p)

    mutation_rate = 0.1
    elite_size = max(1, min(5, population // 4))
    index_pool = list(range(n))

    pop = [rng.sample(index_pool, n) for _ in range(population)]
    score_cache: dict[tuple[int, ...], dict[str, int]] = {}

    def eval_order(ordering: list[int]) -> tuple[dict[str, int], int]:
        key = tuple(ordering)
        if key not in score_cache:
            score_cache[key] = evaluator.score(ordering)
        scores = score_cache[key]
        return scores, _fitness(scores, objective)

    evaluated = [eval_order(ind) for ind in pop]
    evals = len(score_cache)

    best_order = None
    best_scores = None
    best_fit = None
    history: list[int] = []

    for _ in range(generations):
        fitnesses = [fit for _, fit in evaluated]
        ranked_idxs = sorted(range(len(pop)), key=lambda i: fitnesses[i])
        pop = [pop[i] for i in ranked_idxs]
        evaluated = [evaluated[i] for i in ranked_idxs]
        fitnesses = [fit for _, fit in evaluated]

        if best_fit is None or fitnesses[0] < best_fit:
            best_fit = fitnesses[0]
            best_order = list(pop[0])
            best_scores = dict(evaluated[0][0])
        history.append(fitnesses[0])

        next_pop = [list(ind) for ind in pop[:elite_size]]
        while len(next_pop) < population:
            p1 = _tournament_select(pop, fitnesses, rng, k=3)
            p2 = _tournament_select(pop, fitnesses, rng, k=3)
            c1, c2 = _pmx_crossover(p1, p2, rng)
            c1 = _mutate_swap(c1, rng, mutation_rate)
            c2 = _mutate_swap(c2, rng, mutation_rate)
            next_pop.append(c1)
            if len(next_pop) < population:
                next_pop.append(c2)

        pop = next_pop
        evaluated = [eval_order(ind) for ind in pop]
        evals = len(score_cache)

    # The children of the last generation are evaluated above; keep their best too
    # (before this check the loop discarded them).
    for ind, (scores, fit) in zip(pop, evaluated):
        if best_fit is None or fit < best_fit:
            best_fit = fit
            best_order = list(ind)
            best_scores = dict(scores)

    return {
        "method": "ga_pmx_tournament_elite",
        "objective": objective,
        "evaluations": evals,
        "best_ordering": best_order,
        "best_fitness": best_fit,
        "best_scores": best_scores,
        "history": history,
        "seed": seed,
    }


def exhaustive_search(
    transitions_csv: str,
    n: int,
    p: int,
    objective: str = "minimize_mu_E",
    max_evaluations: int | None = None,
) -> dict[str, Any]:
    evaluator = _ReferenceEvaluator(transitions_csv, n, p)

    best_order = None
    best_fit = None
    best_scores = None
    evals = 0

    for perm in itertools.permutations(range(n)):
        if max_evaluations is not None and evals >= max_evaluations:
            break
        scores = evaluator.score(list(perm))
        fit = _fitness(scores, objective)
        evals += 1
        if best_fit is None or fit < best_fit:
            best_fit = fit
            best_order = list(perm)
            best_scores = scores

    return {
        "method": "native_exhaustive",
        "objective": objective,
        "evaluations": evals,
        "best_ordering": best_order,
        "best_fitness": best_fit,
        "best_scores": best_scores,
        "truncated": max_evaluations is not None and evals >= max_evaluations,
    }


def branch_and_bound_search(
    transitions_csv: str,
    n: int,
    p: int,
    objective: str = "minimize_mu_E",
    max_evaluations: int | None = None,
) -> dict[str, Any]:
    evaluator = _ReferenceEvaluator(transitions_csv, n, p)
    used = [False] * n
    current = [0] * n
    evals = 0
    pruned = 0

    best_fit = float("inf")
    best_order: list[int] | None = None
    best_scores: dict[str, int] | None = None

    # Cota del objetivo (fitness = mu_E o -mu_A) válida para toda completación del prefijo.
    if objective == "minimize_mu_E":
        bound = evaluator.partial_mu_e
    else:
        bound = lambda prefix: -evaluator.mu_a_upper_bound(prefix)  # noqa: E731

    def backtrack(depth: int) -> None:
        nonlocal evals, pruned, best_fit, best_order, best_scores

        if max_evaluations is not None and evals >= max_evaluations:
            return

        if depth == n:
            scores = evaluator.score(current)
            fit = _fitness(scores, objective)
            evals += 1
            if fit < best_fit:
                best_fit = fit
                best_order = list(current)
                best_scores = scores
            return

        for g in range(n):
            if used[g]:
                continue
            current[depth] = g
            if bound(current[: depth + 1]) > best_fit:
                pruned += 1
                continue
            used[g] = True
            backtrack(depth + 1)
            used[g] = False

    backtrack(0)
    return {
        "method": "branch_and_bound_partial_muE" if objective == "minimize_mu_E" else "branch_and_bound_upper_muA",
        "objective": objective,
        "evaluations": evals,
        "pruned": pruned,
        "best_ordering": best_order,
        "best_fitness": best_fit if best_fit != float("inf") else None,
        "best_scores": best_scores,
        "truncated": max_evaluations is not None and evals >= max_evaluations,
    }
