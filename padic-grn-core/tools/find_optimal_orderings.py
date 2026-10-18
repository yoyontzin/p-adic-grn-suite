"""
Búsqueda exhaustiva de ordenamientos óptimos para C. elegans y Drosophila.
Algoritmo little-endian, idéntico a padic_grn_core/compute_mu.py: niveles
n = 0, ..., N-1, con mu_E + mu_A + mu_I = N * p^N.

Uso:
    python3 find_optimal_orderings.py
"""
from itertools import permutations


def config_to_vector(config, N, p):
    """Little-endian: v[0]=LSB."""
    v = []
    t = config
    for _ in range(N):
        v.append(t % p)
        t //= p
    return v


def vector_to_config(v, p):
    return sum(v[i] * (p ** i) for i in range(len(v)))


def compute_mu(tt, N, p, A):
    """Calcula (mu_E, mu_A, mu_I) para el ordenamiento A (niveles n = 0..N-1)."""
    mu_E = mu_A = mu_I = 0
    for n in range(0, N):
        nb = p ** n
        sz = p ** (N - n)
        w  = sz  # peso de cada bola en nivel n
        for bi in range(nb):
            imgs = []
            for k in range(sz):
                y = [0] * N
                e = bi
                for j in range(n):
                    y[j] = e % p; e //= p
                e = k
                for j in range(n, N):
                    y[j] = e % p; e //= p
                x = [0] * N
                for j in range(N):
                    x[A[j]] = y[j]
                cfg = vector_to_config(x, p)
                img_cfg = tt[cfg]
                img_x = config_to_vector(img_cfg, N, p)
                imgs.append([img_x[A[j]] for j in range(N)])
            # primer coordenada que diverge
            M = N
            for j in range(N):
                if any(imgs[k][j] != imgs[0][j] for k in range(1, len(imgs))):
                    M = j
                    break
            if   M < n:  mu_E += w
            elif M == n: mu_I += w
            else:        mu_A += w
    return mu_E, mu_A, mu_I


# ──────────────────────────────────────────────
# Tablas de transición (de padic-grn-demo.html)
# ──────────────────────────────────────────────

TT_DROSO = [
    0,3,0,7,0,1,0,5,0,3,0,3,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,
    56,59,56,63,32,33,32,37,40,43,40,43,32,33,32,33,
    48,49,48,49,32,33,32,33,32,33,32,33,32,33,32,33
]

TT_CELEGANS = [
    155,219,23,87,155,219,23,87,186,250,54,118,186,250,54,118,
    187,251,55,119,187,251,55,119,154,218,22,86,154,218,22,86,
    155,219,23,87,155,219,23,87,186,250,54,118,186,250,54,118,
    187,251,55,119,187,251,55,119,154,218,22,86,154,218,22,86,
    43,107,39,103,43,107,39,103,10,74,6,70,10,74,6,70,
    11,75,7,71,11,75,7,71,42,106,38,102,42,106,38,102,
    43,107,39,103,43,107,39,103,10,74,6,70,10,74,6,70,
    11,75,7,71,11,75,7,71,42,106,38,102,42,106,38,102,
    152,152,20,20,152,152,20,20,185,185,53,53,185,185,53,53,
    184,184,52,52,184,184,52,52,153,153,21,21,153,153,21,21,
    152,152,20,20,152,152,20,20,185,185,53,53,185,185,53,53,
    184,184,52,52,184,184,52,52,153,153,21,21,153,153,21,21,
    40,40,36,36,40,40,36,36,9,9,5,5,9,9,5,5,
    8,8,4,4,8,8,4,4,41,41,37,37,41,41,37,37,
    40,40,36,36,40,40,36,36,9,9,5,5,9,9,5,5,
    8,8,4,4,8,8,4,4,41,41,37,37,41,41,37,37
]


def search(name, tt, N, p, known_perm=None):
    total = 1
    for i in range(1, N + 1):
        total *= i
    print(f"\n{'='*60}")
    print(f"{name}  N={N}, p={p}  ({total} permutaciones)")
    print(f"{'='*60}")

    if known_perm:
        muE, muA, muI = compute_mu(tt, N, p, known_perm)
        print(f"Valor actual en app  A={known_perm}")
        print(f"  mu_E={muE}  mu_A={muA}  mu_I={muI}")

    best_E = float('inf')
    best_perm = None
    best_A_val = 0
    best_I_val = 0

    for i, perm in enumerate(permutations(range(N))):
        if i % 5000 == 0:
            print(f"  [{i}/{total}] mejor hasta ahora: mu_E={best_E}", flush=True)
        A = list(perm)
        muE, muA, muI = compute_mu(tt, N, p, A)
        if muE < best_E:
            best_E = muE
            best_A_val = muA
            best_I_val = muI
            best_perm = A[:]

    inv = N * (p ** N)
    print(f"\n  ── Óptimo encontrado ──")
    print(f"  A     = {best_perm}")
    print(f"  mu_E  = {best_E}")
    print(f"  mu_A  = {best_A_val}")
    print(f"  mu_I  = {best_I_val}")
    print(f"  Σ = {best_E+best_A_val+best_I_val}  (debe ser {inv})")


if __name__ == "__main__":
    search(
        "Drosophila melanogaster",
        TT_DROSO, N=6, p=2,
        known_perm=[0, 1, 2, 3, 4, 5]
    )
    search(
        "C. elegans",
        TT_CELEGANS, N=8, p=2,
        known_perm=[1, 0, 2, 3, 6, 5, 4, 7]
    )
