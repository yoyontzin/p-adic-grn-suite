"""CLI funcional nativa para padic-grn-search."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .engine import branch_and_bound_search, random_ga_search


def _write_json(path: str, payload: dict) -> Path:
    out = Path(path).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
    return out


def main_ga() -> None:
    parser = argparse.ArgumentParser(description="Búsqueda GA (torneo + PMX + mutación + élite).")
    parser.add_argument("transitions_csv")
    parser.add_argument("N", type=int)
    parser.add_argument("p", type=int)
    parser.add_argument("--objective", choices=["minimize_mu_E", "maximize_mu_A"], default="minimize_mu_E")
    parser.add_argument("--generations", type=int, default=50)
    parser.add_argument("--population", type=int, default=30)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--output-json", default="search_ga_result.json")
    args = parser.parse_args()
    try:
        result = random_ga_search(
            transitions_csv=args.transitions_csv,
            n=args.N,
            p=args.p,
            objective=args.objective,
            generations=args.generations,
            population=args.population,
            seed=args.seed,
        )
        out = _write_json(args.output_json, result)
        print(f"OK: resultado GA guardado en {out}")
    except Exception as exc:  # pragma: no cover - CLI guard
        print(f"ERROR pgsearch-ga: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


def main_bb_min_mue() -> None:
    parser = argparse.ArgumentParser(description="Búsqueda Branch-and-Bound para minimizar mu_E.")
    parser.add_argument("transitions_csv")
    parser.add_argument("N", type=int)
    parser.add_argument("p", type=int)
    parser.add_argument("--max-evaluations", type=int, default=None)
    parser.add_argument("--output-json", default="search_bb_min_result.json")
    args = parser.parse_args()
    try:
        result = branch_and_bound_search(
            transitions_csv=args.transitions_csv,
            n=args.N,
            p=args.p,
            objective="minimize_mu_E",
            max_evaluations=args.max_evaluations,
        )
        out = _write_json(args.output_json, result)
        print(f"OK: resultado BB-min guardado en {out}")
    except Exception as exc:  # pragma: no cover - CLI guard
        print(f"ERROR pgsearch-bb-min-muE: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


def main_bb_max_mua() -> None:
    parser = argparse.ArgumentParser(description="Búsqueda Branch-and-Bound para maximizar mu_A.")
    parser.add_argument("transitions_csv")
    parser.add_argument("N", type=int)
    parser.add_argument("p", type=int)
    parser.add_argument("--max-evaluations", type=int, default=None)
    parser.add_argument("--output-json", default="search_bb_max_result.json")
    args = parser.parse_args()
    try:
        result = branch_and_bound_search(
            transitions_csv=args.transitions_csv,
            n=args.N,
            p=args.p,
            objective="maximize_mu_A",
            max_evaluations=args.max_evaluations,
        )
        out = _write_json(args.output_json, result)
        print(f"OK: resultado BB-max guardado en {out}")
    except Exception as exc:  # pragma: no cover - CLI guard
        print(f"ERROR pgsearch-bb-max-muA: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
