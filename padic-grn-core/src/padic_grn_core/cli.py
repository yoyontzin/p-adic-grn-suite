"""CLI funcional para padic-grn-core."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any


def _workspace_root() -> Path:
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "20250115-p-adic-GRN-rewrite").exists():
            return parent
    raise RuntimeError("No se encontró la raíz del workspace.")


def _load_core_module() -> Any:
    root = _workspace_root()
    module_path = (
        root
        / "20250115-p-adic-GRN-rewrite"
        / "articulo-final"
        / "supplementary"
        / "code"
        / "core"
        / "compute_mu_generalized.py"
    )
    if not module_path.exists():
        raise FileNotFoundError(f"No existe el módulo core: {module_path}")
    spec = importlib.util.spec_from_file_location("core_compute_mu_generalized", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("No se pudo cargar el módulo core de referencia.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _parse_ordering(ordering_text: str, n: int) -> list[int]:
    order = [int(x) for x in ordering_text.replace(",", " ").split()]
    if sorted(order) != list(range(n)):
        raise ValueError(f"ordering inválido; debe ser permutación de 0..{n-1}")
    return order


def _write_json(output_path: Path, payload: dict[str, Any]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def main_score_ordering() -> None:
    parser = argparse.ArgumentParser(
        description="Calcula mu_E, mu_A, mu_I para un ordering dado."
    )
    parser.add_argument("transitions_csv", help="CSV con columnas config_from,config_to")
    parser.add_argument("N", type=int, help="Número de variables")
    parser.add_argument("p", type=int, help="Base discreta")
    parser.add_argument("ordering", help='Permutación, ej. "0 1 2 3" o "0,1,2,3"')
    parser.add_argument(
        "--output-json",
        default="scores.json",
        help="Archivo de salida JSON (default: scores.json)",
    )
    args = parser.parse_args()

    try:
        module = _load_core_module()
        table = module.read_transition_table_csv(args.transitions_csv, args.p)
        ordering = _parse_ordering(args.ordering, args.N)
        scores = module.compute_scores(args.p, args.N, ordering, table, return_profiles=True)
        payload = {
            "software": "padic-grn-core",
            "version": "0.1.0",
            "N": args.N,
            "p": args.p,
            "ordering": ordering,
            "scores": scores,
        }
        out = Path(args.output_json).resolve()
        _write_json(out, payload)
        print(f"OK: resultados guardados en {out}")
    except Exception as exc:  # pragma: no cover - CLI guard
        print(f"ERROR pgcore-score-ordering: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


def main_batch_scores() -> None:
    parser = argparse.ArgumentParser(
        description="Calcula scores para una lista de orderings (batch)."
    )
    parser.add_argument("transitions_csv", help="CSV con transiciones")
    parser.add_argument("N", type=int, help="Número de variables")
    parser.add_argument("p", type=int, help="Base discreta")
    parser.add_argument(
        "orderings_csv",
        help="CSV con columna ordering (string con enteros separados por espacio/coma)",
    )
    parser.add_argument(
        "--output-csv",
        default="batch_scores.csv",
        help="Archivo de salida CSV (default: batch_scores.csv)",
    )
    args = parser.parse_args()

    try:
        module = _load_core_module()
        table = module.read_transition_table_csv(args.transitions_csv, args.p)

        rows_out: list[dict[str, Any]] = []
        with open(args.orderings_csv, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if "ordering" not in (reader.fieldnames or []):
                raise ValueError("orderings_csv debe contener columna 'ordering'")
            for row in reader:
                ordering_raw = row["ordering"].strip()
                ordering = _parse_ordering(ordering_raw, args.N)
                scores = module.compute_scores(args.p, args.N, ordering, table, return_profiles=False)
                rows_out.append(
                    {
                        "ordering": " ".join(str(x) for x in ordering),
                        "mu_E": scores["mu_E"],
                        "mu_A": scores["mu_A"],
                        "mu_I": scores["mu_I"],
                    }
                )

        out_path = Path(args.output_csv).resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["ordering", "mu_E", "mu_A", "mu_I"])
            writer.writeheader()
            writer.writerows(rows_out)
        print(f"OK: batch guardado en {out_path} ({len(rows_out)} filas)")
    except Exception as exc:  # pragma: no cover - CLI guard
        print(f"ERROR pgcore-batch-scores: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
