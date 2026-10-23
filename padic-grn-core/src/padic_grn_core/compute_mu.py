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
from typing import Dict, List, Set, Tuple, Optional, Sequence
from functools import lru_cache
import multiprocessing as mp
from collections import defaultdict, OrderedDict

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
# EXACT FAST PATH FOR p = 2 (BIT FOLDING)
# ============================================================================
# Exact evaluation of the scores mu_E, mu_A, mu_I for Boolean networks (p = 2),
# with bit operations and standard-library Python only.
#
# It returns exactly the same integers as compute_scores of the reference engine
# (compute_mu_generalized.py) for every ordering, provided the transition table is
# complete (2^N entries, configurations and images in 0..2^N-1).
#
# Ordering and balls (little-endian, as in the reference engine). Configuration c
# encodes gene g in bit g. An ordering A = (a_0, ..., a_{N-1}) lists the genes from
# the most significant coordinate. A ball of level n fixes the genes a_0, ..., a_{n-1}
# and leaves the other N-n genes free. Let OR and AND be the bitwise OR and AND of the
# images F(x) over the ball, and V = OR ^ AND the genes whose next state is not
# constant on the ball. The reference engine takes M = the first position j of the
# ordering with gene a_j in V and classifies the ball as E (M < n), I (M = n) or
# A (M > n). Hence
#     E  <=>  V meets {a_0, ..., a_{n-1}},
#     I  <=>  V misses {a_0, ..., a_{n-1}} and contains a_n,
#     A  otherwise,
# and the level-n counts depend only on the set {a_0, ..., a_{n-1}} and on a_n.
#
# Algorithm. Relabel the configurations so that bit j of the index y holds gene a_j
# (the images keep the original labels). Then a ball of level n is a residue class of
# y modulo 2^n, and the OR/AND over the free genes is obtained by folding the top bit
# of the array, one level at a time, from n = N-1 down to n = 0. Total cost O(2^N)
# list operations per ordering, against O(N^2 2^N) interpreted operations in the
# reference engine.
#
# Level conventions. first_level = 0 sums n = 0, ..., N-1 (mu_E + mu_A + mu_I = N 2^N,
# current convention). first_level = 1 sums n = 1, ..., N-1 (legacy convention of the
# 2026-04 registration, mu_E + mu_A + mu_I = (N-1) 2^N). mu_E does not depend on the
# convention, because the root ball (n = 0) is never expanding. The root ball is A when
# the next state of a_0 is constant on the whole space and I otherwise, so the legacy
# convention removes its weight 2^N from mu_I, or from mu_A in that constant case.
def table_supports_fast_path(N: int, p: int, transition_table: Dict[int, int]) -> bool:
    """True if p = 2 and the table is complete with all values in range."""
    if p != 2 or N < 1 or N > 26:
        return False
    size = 1 << N
    if len(transition_table) != size:
        return False
    try:
        for c in range(size):
            v = transition_table[c]
            if not (0 <= v < size):
                return False
    except (KeyError, TypeError):
        return False
    return True


class FastMuP2:
    """Exact scores for p = 2 by bit folding.

    scores(A)            -> dict(mu_E, mu_A, mu_I[, E, A, I])
    swap_scores(A, i, j) -> scores of A with positions i and j exchanged, computed
                            incrementally from A (levels outside i..j are reused).
    """

    def __init__(self, N: int, transition_table: Dict[int, int], first_level: int = 0,
                 memo_size: int = 1 << 16):
        if first_level not in (0, 1):
            raise ValueError("first_level must be 0 or 1")
        self.N = N
        self.size = 1 << N
        self.F = [int(transition_table[c]) for c in range(self.size)]
        self.first_level = first_level
        self.memo_size = memo_size
        self._memo: "OrderedDict[Tuple[int, ...], Tuple[Tuple[int, ...], Tuple[int, ...], Tuple[int, ...]]]" = OrderedDict()
        # Checkpoint for incremental swaps: ordering and the folded arrays of every level.
        self._ck_order: Optional[Tuple[int, ...]] = None
        self._ck_and: List[List[int]] = []
        self._ck_or: List[List[int]] = []
        self._ck_counts: Optional[Tuple[List[int], List[int], List[int]]] = None
        self.n_full = 0
        self.n_incremental = 0
        self.n_memo_hits = 0

    # ------------------------------------------------------------------ core
    def _relabel(self, A: Sequence[int]) -> List[int]:
        X = [0]
        for g in A:
            b = 1 << g
            X += [x | b for x in X]
        F = self.F
        return [F[x] for x in X]

    @staticmethod
    def _count(AND: List[int], OR: List[int], pmask: int, gbit: int) -> Tuple[int, int]:
        e = i = 0
        for a, o in zip(AND, OR):
            v = a ^ o
            if v & pmask:
                e += 1
            elif v & gbit:
                i += 1
        return e, i

    def _full(self, A: Sequence[int], keep: bool):
        N = self.N
        G = self._relabel(A)
        AND = G
        OR = G
        prefix = [0] * (N + 1)
        for k in range(N):
            prefix[k + 1] = prefix[k] | (1 << A[k])
        E = [0] * N
        I = [0] * N
        Acnt = [0] * N
        ands: List[List[int]] = [None] * (N + 1)  # type: ignore[list-item]
        ors: List[List[int]] = [None] * (N + 1)  # type: ignore[list-item]
        if keep:
            ands[N] = G
            ors[N] = G
        for n in range(N - 1, -1, -1):
            half = 1 << n
            AND = [a & b for a, b in zip(AND[:half], AND[half:])]
            OR = [a | b for a, b in zip(OR[:half], OR[half:])]
            e, i = self._count(AND, OR, prefix[n], 1 << A[n])
            E[n] = e
            I[n] = i
            Acnt[n] = half - e - i
            if keep:
                ands[n] = AND
                ors[n] = OR
        self.n_full += 1
        if keep:
            self._ck_order = tuple(A)
            self._ck_and = ands
            self._ck_or = ors
            self._ck_counts = (E, Acnt, I)
        return E, Acnt, I

    def level_counts(self, A: Sequence[int]) -> Tuple[Tuple[int, ...], Tuple[int, ...], Tuple[int, ...]]:
        """(E_n, A_n, I_n) for n = 0..N-1, memoized by ordering."""
        key = tuple(A)
        hit = self._memo.get(key)
        if hit is not None:
            self._memo.move_to_end(key)
            self.n_memo_hits += 1
            return hit
        E, Acnt, I = self._full(key, keep=False)
        res = (tuple(E), tuple(Acnt), tuple(I))
        self._remember(key, res)
        return res

    def _remember(self, key, res) -> None:
        if self.memo_size <= 0:
            return
        self._memo[key] = res
        if len(self._memo) > self.memo_size:
            self._memo.popitem(last=False)

    def _to_scores(self, counts, return_profiles: bool) -> Dict:
        E, Acnt, I = counts
        N = self.N
        muE = muA = muI = 0
        for n in range(self.first_level, N):
            w = 1 << (N - n)
            muE += E[n] * w
            muA += Acnt[n] * w
            muI += I[n] * w
        out = {"mu_E": muE, "mu_A": muA, "mu_I": muI, "p": 2, "N": N}
        if return_profiles:
            lv = range(self.first_level, N)
            out["E"] = {n: E[n] for n in lv}
            out["A"] = {n: Acnt[n] for n in lv}
            out["I"] = {n: I[n] for n in lv}
        return out

    def scores(self, A: Sequence[int], return_profiles: bool = False) -> Dict:
        return self._to_scores(self.level_counts(A), return_profiles)

    # ------------------------------------------------------- incremental swap
    @staticmethod
    def _swap_index_bits(arr: List[int], i: int, j: int) -> List[int]:
        """New array b with b[y] = arr[y'] where y' is y with bits i < j exchanged.
        len(arr) = 2^(j+1), so bit j is the top bit of the index."""
        h = 1 << j
        L = arr[:h]
        U = arr[h:]
        bs = 1 << i
        newL: List[int] = []
        newU: List[int] = []
        for s in range(0, h, bs << 1):
            # block [s, s+bs) has bit i = 0, block [s+bs, s+2bs) has bit i = 1
            newL += L[s:s + bs]
            newL += U[s:s + bs]
            newU += L[s + bs:s + 2 * bs]
            newU += U[s + bs:s + 2 * bs]
        return newL + newU

    def swap_level_counts(self, A: Sequence[int], i: int, j: int):
        """Level counts of A with positions i and j exchanged (exact)."""
        if i == j:
            return self.level_counts(A)
        if i > j:
            i, j = j, i
        B = list(A)
        B[i], B[j] = B[j], B[i]
        key = tuple(B)
        hit = self._memo.get(key)
        if hit is not None:
            self._memo.move_to_end(key)
            self.n_memo_hits += 1
            return hit
        if self._ck_order != tuple(A):
            self._full(A, keep=True)
        E0, A0, I0 = self._ck_counts  # type: ignore[misc]
        E = list(E0)
        Acnt = list(A0)
        I = list(I0)
        AND = self._swap_index_bits(self._ck_and[j + 1], i, j)
        OR = self._swap_index_bits(self._ck_or[j + 1], i, j)
        pmask = 0
        for k in range(i):
            pmask |= 1 << B[k]
        prefix = [0] * (j + 1)
        acc = pmask
        for n in range(i, j + 1):
            prefix[n] = acc
            acc |= 1 << B[n]
        for n in range(j, i - 1, -1):
            half = 1 << n
            AND = [a & b for a, b in zip(AND[:half], AND[half:])]
            OR = [a | b for a, b in zip(OR[:half], OR[half:])]
            e, ii = self._count(AND, OR, prefix[n], 1 << B[n])
            E[n] = e
            I[n] = ii
            Acnt[n] = half - e - ii
        self.n_incremental += 1
        res = (tuple(E), tuple(Acnt), tuple(I))
        self._remember(key, res)
        return res

    def swap_scores(self, A: Sequence[int], i: int, j: int, return_profiles: bool = False) -> Dict:
        return self._to_scores(self.swap_level_counts(A, i, j), return_profiles)

    def neighbor_level_counts(self, A: Sequence[int], B: Sequence[int]):
        """Level counts of an ordering B that differs from A only inside one block of
        positions lo..hi (a transposition, an insertion move, a segment reversal...).
        The levels outside lo..hi are reused from A (exact)."""
        N = self.N
        lo = 0
        while lo < N and A[lo] == B[lo]:
            lo += 1
        if lo == N:
            return self.level_counts(A)
        hi = N - 1
        while A[hi] == B[hi]:
            hi -= 1
        key = tuple(B)
        hit = self._memo.get(key)
        if hit is not None:
            self._memo.move_to_end(key)
            self.n_memo_hits += 1
            return hit
        if self._ck_order != tuple(A):
            self._full(A, keep=True)
        pos = {g: k for k, g in enumerate(A)}
        # index map of level hi+1: bit k of the new index holds gene B[k]
        idx = [0]
        for k in range(hi + 1):
            b = 1 << pos[B[k]]
            idx += [x | b for x in idx]
        oa = self._ck_and[hi + 1]
        oo = self._ck_or[hi + 1]
        AND = [oa[t] for t in idx]
        OR = [oo[t] for t in idx]
        E0, A0, I0 = self._ck_counts  # type: ignore[misc]
        E = list(E0)
        Acnt = list(A0)
        I = list(I0)
        prefix = [0] * (hi + 1)
        acc = 0
        for n in range(hi + 1):
            prefix[n] = acc
            acc |= 1 << B[n]
        for n in range(hi, lo - 1, -1):
            half = 1 << n
            AND = [a & b for a, b in zip(AND[:half], AND[half:])]
            OR = [a | b for a, b in zip(OR[:half], OR[half:])]
            e, ii = self._count(AND, OR, prefix[n], 1 << B[n])
            E[n] = e
            I[n] = ii
            Acnt[n] = half - e - ii
        self.n_incremental += 1
        res = (tuple(E), tuple(Acnt), tuple(I))
        self._remember(key, res)
        return res

    def neighbor_scores(self, A: Sequence[int], B: Sequence[int], return_profiles: bool = False) -> Dict:
        return self._to_scores(self.neighbor_level_counts(A, B), return_profiles)

    def clear(self) -> None:
        self._memo.clear()
        self._ck_order = None
        self._ck_and = []
        self._ck_or = []
        self._ck_counts = None

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
        # Exact fast path for p = 2 with a complete table (same integers as the
        # per-ball computation below, verified against it and against the exact DP).
        self._fast: Optional[FastMuP2] = None
        if table_supports_fast_path(N, p, transition_table):
            self._fast = FastMuP2(N, transition_table, first_level=0)

    def compute_ball_diameter_cached(self, n: int, ball_index: int, A: List[int]) -> int:
        """Compute M (the level such that t_{n,m} = 1/p^M), with cache."""
        cache_key = (n, ball_index, tuple(A))
        if cache_key in self._cache:
            return self._cache[cache_key]
        if len(self._cache) >= 4_000_000:
            # one entry per ball and per ordering: without a bound it grows by
            # about 2^N entries for every new ordering
            self._cache.clear()
        M = compute_ball_diameter_generalized(
            n, ball_index, self.N, self.p, A, self.transition_table
        )
        self._cache[cache_key] = M
        return M

    def compute_mu(self, A: List[int], use_parallel: bool = False) -> int:
        """Compute mu_E (expansion measure)."""
        if self._fast is not None:
            return self._fast.scores(A)["mu_E"]
        return compute_mu_generalized(
            self.N, self.p, A, self.transition_table,
            use_parallel=use_parallel
        )

    def compute_scores_cached(
        self, A: List[int], return_profiles: bool = True
    ) -> Dict:
        """Compute mu_E, mu_A, mu_I and the E/A/I profiles using the cache of M."""
        if self._fast is not None:
            return self._fast.scores(A, return_profiles=return_profiles)
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

    def compute_scores_swap(
        self, A: List[int], i: int, j: int, return_profiles: bool = False
    ) -> Dict:
        """Scores of A with positions i and j exchanged. On the p = 2 fast path the
        levels outside i..j are reused from A, otherwise the swapped ordering is
        evaluated in full. The result equals compute_scores_cached(swapped)."""
        if self._fast is not None:
            return self._fast.swap_scores(A, i, j, return_profiles=return_profiles)
        B = list(A)
        B[i], B[j] = B[j], B[i]
        return self.compute_scores_cached(B, return_profiles=return_profiles)

    def compute_scores_neighbor(
        self, A: List[int], B: List[int], return_profiles: bool = False
    ) -> Dict:
        """Scores of an ordering B that differs from A in one block of positions
        (insertion move, transposition, reversal). On the p = 2 fast path the
        levels outside that block are reused from A. Equals compute_scores_cached(B)."""
        if self._fast is not None:
            return self._fast.neighbor_scores(A, B, return_profiles=return_profiles)
        return self.compute_scores_cached(list(B), return_profiles=return_profiles)

    def clear_cache(self) -> None:
        self._cache.clear()
        if self._fast is not None:
            self._fast.clear()

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
