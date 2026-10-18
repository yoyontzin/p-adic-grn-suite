#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generalized computation of the stability measure mu for gene regulatory networks.

Vendored in padic-grn-core from the Paper II reproducibility bundle
(supplementary/code/core/compute_mu_generalized.py, commit 681aeee), with the
partial bounds partial_level_counts, partial_mu_A and mu_A_upper_bound added
for branch and bound. Levels are n = 0, ..., N-1 and mu_E + mu_A + mu_I = N * p^N.

This version works for any base p (not only p=2), so that each gene can take
p different states (0, 1, ..., p-1).

Features:
- Supports any p (p=2 for Boolean networks, p=3, p=4, etc.)
- Uses a transition table (more general than polynomial update functions)
- Implements the observable Z_p-radius (rad_Zp) of the image of each cylinder
  through the "first digit of disagreement" between images. This agrees with
  t_{n,m} and lambda_{n,m} as defined in the companion paper (observable
  Z_p-radius).
- Uses memoization and efficient data structures
- Can run in parallel for large networks

Usage:
    python compute_mu.py transitions.csv N p [ordering]

Example (Boolean, p=2), run from the padic-grn-core/ folder:
    python src/padic_grn_core/compute_mu.py examples/transitions_athaliana_N13_p2.csv 13 2 "12 3 4 8 5 1 7 2 10 11 9 6 0"

Example (p=3, each gene has 3 states):
    python src/padic_grn_core/compute_mu.py <transitions_p3.csv> 10 3 "0 1 2 3 4 5 6 7 8 9"
"""

import sys
import csv
import math
from typing import Dict, List, Set, Tuple, Optional
from functools import lru_cache
import multiprocessing as mp
from collections import defaultdict

# ============================================================================
# GENERALIZED HELPER FUNCTIONS
# ============================================================================

def config_to_vector(config: int, N: int, p: int) -> List[int]:
    """Convert a configuration index to its base-p digit vector."""
    vector = []
    temp = config
    for i in range(N):
        vector.append(temp % p)
        temp //= p
    return vector

def vector_to_config(vector: List[int], N: int, p: int) -> int:
    """Convert a base-p digit vector to its configuration index."""
    config = 0
    for i in range(N):
        config += vector[i] * (p ** i)
    return config

def read_transition_table_csv(filename: str, p: int) -> Dict[int, int]:
    """
    Read a transition table from a CSV file.

    Args:
        filename: path to the CSV file
        p: base (number of states per gene)

    Returns:
        Dictionary: current_config -> next_config
    """
    transition_table = {}
    
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            next(reader)  # Skip header
            
            for row in reader:
                if len(row) >= 2:
                    config_from = int(row[0])
                    config_to = int(row[1])
                    transition_table[config_from] = config_to
    except FileNotFoundError:
        print(f"Error: could not open file {filename}")
        sys.exit(1)
    except Exception as e:
        print(f"Error while reading the file: {e}")
        sys.exit(1)
    
    return transition_table

# ============================================================================
# IMAGE-BALL DIAMETER (GENERALIZED TO ANY p)
# ============================================================================

def compute_ball_diameter_generalized(n: int, ball_index: int, N: int, p: int, 
                                     A: List[int], transition_table: Dict[int, int]) -> int:
    """
    Compute the diameter of the image ball for B(1/p^n, ball_index).

    Args:
        n: resolution level
        ball_index: index of the ball (0 to p^n - 1)
        N: number of genes
        p: base (number of states per gene)
        A: gene ordering
        transition_table: transition table

    Returns:
        Level M such that t_{n,m} = 1/p^M. Classification: M < n -> E (expanding, Lambda > 1),
        M == n -> I (indifferent, Lambda = 1), M > n -> A (contracting, Lambda < 1).
    """
    # Collect full p-adic digit vectors (length N) for all images in B_{n,m}
    image_vectors: List[List[int]] = []
    num_configs = p ** (N - n)

    for k in range(num_configs):
        y = [0] * N
        e = ball_index
        for j in range(n):
            y[j] = e % p
            e //= p
        e = k
        for j in range(n, N):
            y[j] = e % p
            e //= p
        x = [0] * N
        for j in range(N):
            x[A[j]] = y[j]
        config = vector_to_config(x, N, p)
        if config not in transition_table:
            continue
        next_config = transition_table[config]
        next_x = config_to_vector(next_config, N, p)
        y_image = [next_x[A[j]] for j in range(N)]
        image_vectors.append(y_image)

    if not image_vectors:
        return N

    # M = smallest j in {0,...,N-1} such that two images differ at digit j. M = N if all agree
    M = N
    for j in range(N):
        digits_j = {vec[j] for vec in image_vectors}
        if len(digits_j) > 1:
            M = j
            break
    return M


def classify_ball_with_prefix(
    n: int,
    ball_index: int,
    prefix: List[int],
    N: int,
    p: int,
    transition_table: Dict[int, int],
) -> str:
    """
    Classify the ball B(1/p^n, ball_index) as 'E', 'I' or 'A' using only the
    prefix A[0],...,A[n] of the ordering (len(prefix) >= n+1).

    The ball fixes the genes prefix[0],...,prefix[n-1]. The images are compared
    only on coordinates 0..n of the ordered space (prefix[0],...,prefix[n]), so
    E_n depends only on the prefix of length n+1.

    Returns:
        'E' if M < n, 'I' if M == n, 'A' if M > n.
    """
    assert len(prefix) >= n + 1, "prefix must have at least n+1 entries for level n"
    # Genes fixed by the ball: prefix[0],...,prefix[n-1]
    fixed_set = set(prefix[:n])
    free_list = [g for g in range(N) if g not in fixed_set]
    num_configs = p ** (N - n)
    image_prefix_vectors: List[List[int]] = []  # only coordinates 0..n, in order

    for ki in range(num_configs):
        x = [0] * N
        e = ball_index
        for j in range(n):
            x[prefix[j]] = e % p
            e //= p
        e = ki
        for j in range(len(free_list)):
            x[free_list[j]] = e % p
            e //= p
        config = vector_to_config(x, N, p)
        if config not in transition_table:
            continue
        next_x = config_to_vector(transition_table[config], N, p)
        y_prefix = [next_x[prefix[j]] for j in range(n + 1)]
        image_prefix_vectors.append(y_prefix)

    if not image_prefix_vectors:
        return "A"

    M = n + 1
    for j in range(n + 1):
        digits_j = {vec[j] for vec in image_prefix_vectors}
        if len(digits_j) > 1:
            M = j
            break
    if M < n:
        return "E"
    if M == n:
        return "I"
    return "A"


def partial_level_counts(
    prefix: List[int],
    N: int,
    p: int,
    transition_table: Dict[int, int],
) -> Dict[int, Tuple[int, int, int]]:
    """
    Exact counts (E_n, A_n, I_n) for the levels n = 0, ..., k fixed by the prefix
    (a_0,...,a_k) of an ordering. Level n uses prefix[0..n], so these counts are
    the same for every completion of the prefix.

    Returns:
        dict n -> (E_n, A_n, I_n), with E_n + A_n + I_n = p^n. Empty if the prefix is empty.
    """
    counts: Dict[int, Tuple[int, int, int]] = {}
    for n in range(0, len(prefix)):
        E_n, A_n, I_n = 0, 0, 0
        for i in range(p ** n):
            cls = classify_ball_with_prefix(n, i, prefix, N, p, transition_table)
            if cls == "E":
                E_n += 1
            elif cls == "I":
                I_n += 1
            else:
                A_n += 1
        counts[n] = (E_n, A_n, I_n)
    return counts


def partial_mu_E(
    prefix: List[int],
    N: int,
    p: int,
    transition_table: Dict[int, int],
) -> int:
    """
    Partial sum mu_k = sum_{n=0}^{k} E_n * p^{N-n} for the prefix (a_0,...,a_k).

    prefix has length k+1 (indices 0..k). Level n uses prefix[0..n].
    Bound: mu_E(full ordering) >= partial_mu_E(prefix) for every completion.

    Returns:
        mu_k (int >= 0). If the prefix is empty, returns 0 (no level).
    """
    counts = partial_level_counts(prefix, N, p, transition_table)
    return sum(E_n * p ** (N - n) for n, (E_n, _, _) in counts.items())


def partial_mu_A(
    prefix: List[int],
    N: int,
    p: int,
    transition_table: Dict[int, int],
) -> int:
    """
    Partial sum sum_{n=0}^{k} A_n * p^{N-n} for the prefix (a_0,...,a_k).

    Bound: mu_A(full ordering) >= partial_mu_A(prefix) for every completion.
    """
    counts = partial_level_counts(prefix, N, p, transition_table)
    return sum(A_n * p ** (N - n) for n, (_, A_n, _) in counts.items())


def mu_A_upper_bound(
    prefix: List[int],
    N: int,
    p: int,
    transition_table: Dict[int, int],
) -> int:
    """
    Upper bound on mu_A over every completion of the prefix (a_0,...,a_k):

        mu_A <= partial_mu_A(prefix) + (N - 1 - k) * p^N,

    because each of the N-1-k levels n = k+1, ..., N-1 not yet fixed contributes
    at most p^n * p^{N-n} = p^N. For the empty prefix (k = -1) the bound is N * p^N.
    """
    k = len(prefix) - 1
    return partial_mu_A(prefix, N, p, transition_table) + (N - 1 - k) * p ** N


# ============================================================================
# COMPUTATION OF mu AND THE SCORES mu_E, mu_A, mu_I (GENERALIZED TO ANY p)
# ============================================================================

def compute_scores(
    p: int,
    N: int,
    ordering: List[int],
    transition_table: Dict[int, int],
    return_profiles: bool = True,
    return_lambda_stats: bool = False,
    use_parallel: bool = False,
    num_processes: Optional[int] = None,
) -> Dict:
    """
    Compute the three scores mu_E, mu_A, mu_I and the per-level profiles E(n), A(n), I(n).

    Implements rad_Zp (observable radius) through congruences (first digit of disagreement).
    Classification: M < n -> E (expanding), M = n -> I (indifferent),
    M > n or zero diameter -> A (contracting).

    The levels are n = 0, ..., N-1, so mu_E + mu_A + mu_I = N * p^N.

    Args:
        p: base (number of states per gene).
        N: number of genes.
        ordering: gene ordering (permutation of [0, 1, ..., N-1]).
        transition_table: transition table config -> config.
        return_profiles: if True, include E, A, I per level.
        return_lambda_stats: accepted for interface compatibility. It is not used.
        use_parallel: if True, use parallel processing when there are many cylinders.
        num_processes: number of processes (None = automatic).

    Returns:
        dict with: mu_E, mu_A, mu_I (int), E, A, I (dict n -> count), p, N.
    """
    E_profile: Dict[int, int] = {}
    A_profile: Dict[int, int] = {}
    I_profile: Dict[int, int] = {}
    mu_E, mu_A, mu_I = 0, 0, 0

    for n in range(0, N):  # n = 0, ..., N-1
        E_n, A_n, I_n = 0, 0, 0
        num_balls = p ** n

        if use_parallel and num_balls > 100:
            if num_processes is None:
                num_processes = mp.cpu_count()
            with mp.Pool(processes=num_processes) as pool:
                args = [(n, i, N, p, ordering, transition_table) for i in range(num_balls)]
                results = pool.starmap(compute_ball_diameter_generalized, args)
                for M in results:
                    if M < n:
                        E_n += 1
                    elif M == n:
                        I_n += 1
                    else:
                        A_n += 1
        else:
            for i in range(num_balls):
                M = compute_ball_diameter_generalized(n, i, N, p, ordering, transition_table)
                if M < n:
                    E_n += 1
                elif M == n:
                    I_n += 1
                else:
                    A_n += 1

        assert E_n + A_n + I_n == num_balls, f"E(n)+A(n)+I(n) = p^n at n={n}"
        E_profile[n] = E_n
        A_profile[n] = A_n
        I_profile[n] = I_n
        w = p ** (N - n)
        mu_E += E_n * w
        mu_A += A_n * w
        mu_I += I_n * w

    out: Dict = {
        "mu_E": mu_E,
        "mu_A": mu_A,
        "mu_I": mu_I,
        "p": p,
        "N": N,
    }
    if return_profiles:
        out["E"] = E_profile
        out["A"] = A_profile
        out["I"] = I_profile
    return out


def compute_mu_generalized(N: int, p: int, A: List[int], 
                          transition_table: Dict[int, int],
                          use_parallel: bool = False,
                          num_processes: Optional[int] = None) -> int:
    """
    Compute the stability measure mu_E (expansion) for a given ordering.
    Equivalent to compute_scores(...)['mu_E']. Kept for backward compatibility.
    """
    scores = compute_scores(p, N, A, transition_table, return_profiles=False,
                            use_parallel=use_parallel, num_processes=num_processes)
    return scores["mu_E"]

# ============================================================================
# OPTIMIZATION: MEMOIZATION AND CACHE
# ============================================================================

class MuCalculator:
    """
    Calculator for mu_E (and optionally mu_A, mu_I) with memoization.
    compute_mu returns mu_E (the measure mu of the manuscript).
    """

    def __init__(self, N: int, p: int, transition_table: Dict[int, int]):
        self.N = N
        self.p = p
        self.transition_table = transition_table
        self._cache: Dict[Tuple[int, int, Tuple[int, ...]], int] = {}

    def compute_ball_diameter_cached(self, n: int, ball_index: int, A: List[int]) -> int:
        """Compute M (the level such that t_{n,m} = 1/p^M), with cache."""
        cache_key = (n, ball_index, tuple(A))
        if cache_key in self._cache:
            return self._cache[cache_key]
        M = compute_ball_diameter_generalized(
            n, ball_index, self.N, self.p, A, self.transition_table
        )
        self._cache[cache_key] = M
        return M

    def compute_mu(self, A: List[int], use_parallel: bool = False) -> int:
        """Compute mu_E (expansion measure)."""
        return compute_mu_generalized(
            self.N, self.p, A, self.transition_table,
            use_parallel=use_parallel
        )

    def compute_scores_cached(
        self, A: List[int], return_profiles: bool = True
    ) -> Dict:
        """Compute mu_E, mu_A, mu_I and the E/A/I profiles using the cache of M."""
        E_profile = {}
        A_profile = {}
        I_profile = {}
        mu_E = mu_A = mu_I = 0
        for n in range(0, self.N):
            E_n = A_n = I_n = 0
            for i in range(self.p ** n):
                M = self.compute_ball_diameter_cached(n, i, A)
                if M < n:
                    E_n += 1
                elif M == n:
                    I_n += 1
                else:
                    A_n += 1
            E_profile[n] = E_n
            A_profile[n] = A_n
            I_profile[n] = I_n
            w = self.p ** (self.N - n)
            mu_E += E_n * w
            mu_A += A_n * w
            mu_I += I_n * w
        out = {"mu_E": mu_E, "mu_A": mu_A, "mu_I": mu_I, "p": self.p, "N": self.N}
        if return_profiles:
            out["E"] = E_profile
            out["A"] = A_profile
            out["I"] = I_profile
        return out

    def clear_cache(self) -> None:
        self._cache.clear()

# ============================================================================
# MAIN
# ============================================================================

def main():
    if len(sys.argv) < 4:
        print("Usage: python compute_mu.py <transitions.csv> <N> <p> [ordering] [--parallel]")
        print("\nArguments:")
        print("  transitions.csv: CSV file with the transition table")
        print("  N: number of genes")
        print("  p: base (number of states per gene, p=2 for Boolean networks)")
        print("  ordering: gene ordering (optional)")
        print("  --parallel: use parallel processing")
        print("\nExamples (run from the padic-grn-core/ folder):")
        print("  # Boolean (p=2):")
        print('  python src/padic_grn_core/compute_mu.py examples/transitions_athaliana_N13_p2.csv 13 2 "12 3 4 8 5 1 7 2 10 11 9 6 0"')
        print("\n  # Network with 3 states per gene (p=3):")
        print('  python src/padic_grn_core/compute_mu.py <transitions_p3.csv> 10 3 "0 1 2 3 4 5 6 7 8 9" --parallel')
        sys.exit(1)
    
    filename = sys.argv[1]
    N = int(sys.argv[2])
    p = int(sys.argv[3])
    
    if N <= 0 or N > 30:
        print(f"Error: N must be between 1 and 30 (received: {N})")
        sys.exit(1)
    
    if p < 2 or p > 10:
        print(f"Error: p must be between 2 and 10 (received: {p})")
        sys.exit(1)
    
    use_parallel = '--parallel' in sys.argv
    num_processes = None  # None = automatic if use_parallel

    # Read the transition table
    print(f"Reading transition table from {filename}...")
    transition_table = read_transition_table_csv(filename, p)
    
    expected_configs = p ** N
    if len(transition_table) != expected_configs:
        print(f"Warning: read {len(transition_table)} transitions, expected {expected_configs}")
    
    print(f"Transitions read: {len(transition_table)} of {expected_configs} possible")
    print(f"Base: p = {p} (each gene has {p} states: 0, 1, ..., {p-1})")
    
    # Read the gene ordering
    ordering_arg = None
    for i, arg in enumerate(sys.argv[4:], 4):
        if arg != '--parallel':
            ordering_arg = arg
            break
    
    if ordering_arg:
        A = [int(x) for x in ordering_arg.split()]
    else:
        print(f"\nEnter the gene ordering ({N} integers, separated by spaces):")
        A = [int(x) for x in input().split()]
    
    # Validate the ordering
    if len(A) != N:
        print(f"Error: expected {N} genes, received {len(A)}")
        sys.exit(1)
    
    # Check that it is a valid permutation
    used = [False] * N
    for idx in A:
        if idx < 0 or idx >= N:
            print(f"Error: gene index out of range: {idx}")
            sys.exit(1)
        if used[idx]:
            print(f"Error: duplicate index in the ordering: {idx}")
            sys.exit(1)
        used[idx] = True
    
    # Compute the scores (mu_E, mu_A, mu_I) and the profiles E(n), A(n), I(n)
    print(f"\nComputing scores mu_E, mu_A, mu_I and E/A/I profiles...")
    print(f"Ordering: {A}")
    if use_parallel:
        print(f"Using parallel processing ({mp.cpu_count()} processes)")
    
    scores = compute_scores(p, N, A, transition_table, return_profiles=True,
                            use_parallel=use_parallel, num_processes=num_processes)
    mu_E, mu_A, mu_I = scores["mu_E"], scores["mu_A"], scores["mu_I"]
    
    print(f"\n{'='*50}")
    print(f"  mu_E (expanding)   = {mu_E:,}")
    print(f"  mu_A (contracting) = {mu_A:,}")
    print(f"  mu_I (indifferent) = {mu_I:,}")
    print(f"  mu = mu_E (measure used in the manuscript) = {mu_E:,}")
    print(f"{'='*50}")
    print(f"\nProfiles per level n=0..N-1:")
    print(f"  {'n':<4} {'E(n)':>8} {'A(n)':>8} {'I(n)':>8}  E+A+I=p^n")
    for n in range(0, N):
        E_n, A_n, I_n = scores["E"][n], scores["A"][n], scores["I"][n]
        print(f"  {n:<4} {E_n:>8} {A_n:>8} {I_n:>8}  {E_n + A_n + I_n} = {p**n}")
    print(f"\nParameters:")
    print(f"  N (genes): {N}")
    print(f"  p (states per gene): {p}")
    print(f"  Total number of configurations: {p**N:,}")
    
    return mu_E

if __name__ == "__main__":
    main()
