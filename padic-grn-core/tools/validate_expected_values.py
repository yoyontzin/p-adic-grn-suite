"""Valida los valores de referencia recalculándolos con el motor de padic-grn-core.

Los niveles son n = 0, ..., N-1. El nivel raíz n = 0 tiene una sola bola (todo
el espacio), que nunca es expansiva, de modo que

    mu_E + mu_A + mu_I = N * p^N     y     E(n) + A(n) + I(n) = p^n.

Los valores de A. thaliana y C. elegans coinciden con el bundle de
reproducibilidad de Paper II (mu_dp_exact.js y compute_mu_generalized.py).

Uso (desde padic-grn-core/):
    python3 tools/validate_expected_values.py [--output reporte.md]
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from padic_grn_core import compute_mu  # noqa: E402


@dataclass
class CaseSpec:
    name: str
    csv_path: Path
    ordering: List[int]
    expected_scores: Dict[str, int]
    n: int
    p: int


def check_case(case: CaseSpec) -> Tuple[bool, List[str]]:
    table = compute_mu.read_transition_table_csv(str(case.csv_path), case.p)
    scores = compute_mu.compute_scores(case.p, case.n, case.ordering, table, return_profiles=True)
    issues: List[str] = []

    for key, expected in case.expected_scores.items():
        got = scores.get(key)
        if got != expected:
            issues.append(f"{key}: esperado={expected}, obtenido={got}")

    left = scores["mu_E"] + scores["mu_A"] + scores["mu_I"]
    right = case.n * case.p**case.n
    if left != right:
        issues.append(f"identidad_global: esperado={right}, obtenido={left}")

    for level in range(0, case.n):
        got = scores["E"][level] + scores["A"][level] + scores["I"][level]
        expected = case.p**level
        if got != expected:
            issues.append(f"identidad_nivel_{level}: esperado={expected}, obtenido={got}")

    return (len(issues) == 0), issues


def make_report(cases: List[CaseSpec]) -> Tuple[bool, str]:
    lines = [
        "# Validación de valores de referencia (niveles n = 0..N-1)",
        "",
        "Cada caso se recalcula desde su tabla de transición con padic_grn_core.compute_mu.",
        "",
    ]
    global_ok = True
    for case in cases:
        ok, issues = check_case(case)
        global_ok = global_ok and ok
        lines.append(f"## Caso: {case.name} — {'PASS' if ok else 'FAIL'}")
        lines.append("")
        lines.append(f"- Tabla: `{case.csv_path.relative_to(ROOT)}`")
        lines.append(f"- Ordering: `{case.ordering}`")
        lines.append(f"- Esperados: `{case.expected_scores}`")
        if ok:
            lines.append("- Resultado: todos los checks pasaron.")
        else:
            lines.append("- Issues:")
            lines.extend(f"  - {issue}" for issue in issues)
        lines.append("")
    lines.append(f"## Resultado global: {'PASS' if global_ok else 'FAIL'}")
    lines.append("")
    return global_ok, "\n".join(lines)


def reference_cases() -> List[CaseSpec]:
    ex = ROOT / "examples"
    return [
        CaseSpec(
            name="A. thaliana pi* (N=13, p=2)",
            csv_path=ex / "transitions_athaliana_N13_p2.csv",
            ordering=[12, 3, 4, 8, 5, 1, 7, 2, 10, 11, 9, 6, 0],
            expected_scores={"mu_E": 26776, "mu_A": 53880, "mu_I": 25840},
            n=13,
            p=2,
        ),
        CaseSpec(
            name="C. elegans pi* (N=8, p=2)",
            csv_path=ex / "transitions_celegans_N8_p2.csv",
            ordering=[1, 0, 2, 3, 6, 5, 4, 7],
            expected_scores={"mu_E": 956, "mu_A": 392, "mu_I": 700},
            n=8,
            p=2,
        ),
        CaseSpec(
            name="Drosophila pi* (N=6, p=2)",
            csv_path=ex / "transitions_drosophila_N6_p2.csv",
            ordering=[0, 5, 3, 4, 2, 1],
            expected_scores={"mu_E": 64, "mu_A": 156, "mu_I": 164},
            n=6,
            p=2,
        ),
        CaseSpec(
            name="Toy identidad p=2, N=2",
            csv_path=ex / "transitions_p2_n2_identity.csv",
            ordering=[0, 1],
            expected_scores={"mu_E": 0, "mu_A": 0, "mu_I": 8},
            n=2,
            p=2,
        ),
        CaseSpec(
            name="Toy identidad p=3, N=2",
            csv_path=ex / "transitions_p3_n2_identity.csv",
            ordering=[0, 1],
            expected_scores={"mu_E": 0, "mu_A": 0, "mu_I": 18},
            n=2,
            p=3,
        ),
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description="Valida valores de referencia (n = 0..N-1).")
    parser.add_argument("--output", default=None, help="Archivo Markdown para el reporte (opcional)")
    args = parser.parse_args()

    ok, report = make_report(reference_cases())
    if args.output:
        out = Path(args.output).resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(report, encoding="utf-8")
        print(f"Reporte guardado en {out}")
    else:
        print(report)
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
