"""Motor nativo inicial para visualización en padic-grn-viz."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_scores_json(path: str) -> dict[str, Any]:
    src = Path(path).resolve()
    with src.open("r", encoding="utf-8") as f:
        payload = json.load(f)
    return payload


def plot_score_bars(payload: dict[str, Any], output_path: str) -> Path:
    scores = payload.get("scores", payload.get("best_scores", {}))
    mu_e = scores.get("mu_E", 0)
    mu_a = scores.get("mu_A", 0)
    mu_i = scores.get("mu_I", 0)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(["mu_E", "mu_A", "mu_I"], [mu_e, mu_a, mu_i], color=["#2E86AB", "#6A994E", "#7D53DE"])
    ax.set_title("Resumen de métricas p-ádicas")
    ax.set_ylabel("Valor")
    fig.tight_layout()

    out = Path(output_path).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def plot_ordering_positions(payload: dict[str, Any], output_path: str) -> Path:
    ordering = payload.get("ordering", payload.get("best_ordering", []))
    if not ordering:
        ordering = []

    fig, ax = plt.subplots(figsize=(7, 3))
    if ordering:
        ax.plot(range(len(ordering)), ordering, marker="o")
    ax.set_title("Posición vs índice de gen (ordering)")
    ax.set_xlabel("Posición")
    ax.set_ylabel("Índice de gen")
    fig.tight_layout()

    out = Path(output_path).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out
