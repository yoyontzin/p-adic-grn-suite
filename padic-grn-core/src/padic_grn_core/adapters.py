"""Adaptadores de entrada para padic-grn-core."""

from __future__ import annotations

from typing import Any


def from_grn_table(transitions: dict[int, int], n: int, p: int) -> dict[str, Any]:
    """Adaptador oficial para GRN discretas deterministas.

    Retorna una estructura normalizada para el núcleo de cálculo.
    """
    expected = p**n
    return {
        "kind": "grn_table",
        "N": n,
        "p": p,
        "expected_states": expected,
        "num_transitions": len(transitions),
        "transitions": transitions,
    }


def from_discrete_system(
    transitions: dict[int, int],
    n: int,
    p: int,
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Adaptador para sistemas discretos generales.

    Esta interfaz asume codificación a espacio uniforme p^N.
    """
    payload = from_grn_table(transitions, n, p)
    payload["kind"] = "discrete_system"
    payload["metadata"] = metadata or {}
    return payload
