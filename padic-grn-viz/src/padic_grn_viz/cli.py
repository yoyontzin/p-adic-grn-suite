"""CLI funcional nativa para padic-grn-viz."""

from __future__ import annotations

import argparse
import sys

from .engine import load_scores_json, plot_ordering_positions, plot_score_bars


def main_dependency_graph() -> None:
    parser = argparse.ArgumentParser(
        description="Genera figura de barras para mu_E, mu_A, mu_I."
    )
    parser.add_argument("scores_json", help="Archivo JSON de resultados (core/search)")
    parser.add_argument("--output", default="viz_dependency_scores.png", help="Salida PNG")
    args = parser.parse_args()
    try:
        payload = load_scores_json(args.scores_json)
        out = plot_score_bars(payload, args.output)
        print(f"OK: figura guardada en {out}")
    except Exception as exc:  # pragma: no cover - CLI guard
        print(f"ERROR pgviz-dependency-graph: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


def main_hierarchical_graph() -> None:
    parser = argparse.ArgumentParser(
        description="Genera figura de posiciones del ordering."
    )
    parser.add_argument("scores_json", help="Archivo JSON de resultados (core/search)")
    parser.add_argument("--output", default="viz_hierarchical_ordering.png", help="Salida PNG")
    args = parser.parse_args()
    try:
        payload = load_scores_json(args.scores_json)
        out = plot_ordering_positions(payload, args.output)
        print(f"OK: figura guardada en {out}")
    except Exception as exc:  # pragma: no cover - CLI guard
        print(f"ERROR pgviz-hierarchical-graph: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
