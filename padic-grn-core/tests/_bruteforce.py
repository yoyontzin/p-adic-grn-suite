"""Evaluador de fuerza bruta, escrito directamente desde la definición.

Para el ordenamiento A y el nivel n, la bola de un estado x es la tupla de sus
dígitos en los genes A[0], ..., A[n-1]. M es la primera posición j del
ordenamiento en la que las imágenes de la bola no coinciden en el gen A[j]
(M = N si todas coinciden). E si M < n, I si M = n, A si M > n. Niveles
n = 0, ..., N-1, peso p^(N-n) por bola.

No comparte código con padic_grn_core.compute_mu.
"""

from __future__ import annotations

import itertools
import random
from collections import defaultdict


def digits(x: int, N: int, p: int) -> tuple[int, ...]:
    return tuple((x // p**i) % p for i in range(N))


def level_counts(table: dict[int, int], N: int, p: int, A: list[int]) -> dict[int, tuple[int, int, int]]:
    states = [digits(x, N, p) for x in range(p**N)]
    images = [digits(table[x], N, p) for x in range(p**N)]
    out = {}
    for n in range(N):
        balls: dict[tuple[int, ...], list[int]] = defaultdict(list)
        for x, s in enumerate(states):
            balls[tuple(s[A[j]] for j in range(n))].append(x)
        assert len(balls) == p**n
        e = a = i = 0
        for members in balls.values():
            M = N
            for j in range(N):
                if len({images[x][A[j]] for x in members}) > 1:
                    M = j
                    break
            if M < n:
                e += 1
            elif M == n:
                i += 1
            else:
                a += 1
        out[n] = (e, a, i)
    return out


def scores(table: dict[int, int], N: int, p: int, A: list[int]) -> tuple[int, int, int]:
    mu = [0, 0, 0]
    for n, counts in level_counts(table, N, p, A).items():
        for k in range(3):
            mu[k] += counts[k] * p ** (N - n)
    return tuple(mu)


def optima(table: dict[int, int], N: int, p: int) -> tuple[int, int]:
    """(min mu_E, max mu_A) sobre los N! ordenamientos."""
    all_scores = [scores(table, N, p, list(A)) for A in itertools.permutations(range(N))]
    return min(s[0] for s in all_scores), max(s[1] for s in all_scores)


def random_network(N: int, p: int, seed: int, n_constant: int) -> dict[int, int]:
    """Red aleatoria; los primeros n_constant genes (tras barajar) tienen actualización constante."""
    rng = random.Random(seed)
    genes = list(range(N))
    rng.shuffle(genes)
    constant = {g: rng.randrange(p) for g in genes[:n_constant]}
    table = {}
    for x in range(p**N):
        img = [constant[g] if g in constant else rng.randrange(p) for g in range(N)]
        table[x] = sum(v * p**i for i, v in enumerate(img))
    return table


def constant_genes(table: dict[int, int], N: int, p: int) -> set[int]:
    return {g for g in range(N) if len({digits(y, N, p)[g] for y in table.values()}) == 1}


# (N, p, seed, genes constantes). Todas tienen al menos un gen constante,
# salvo dos de control.
CASES = [
    (3, 2, 1, 1), (4, 2, 2, 1), (4, 2, 3, 2), (5, 2, 4, 1), (5, 2, 5, 0),
    (6, 2, 6, 1), (6, 2, 7, 2),
    (2, 3, 8, 1), (3, 3, 9, 1), (4, 3, 10, 1), (4, 3, 11, 0), (5, 3, 12, 2),
]
