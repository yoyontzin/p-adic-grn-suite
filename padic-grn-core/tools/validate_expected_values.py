"""Valida valores esperados contra evidencia JSON generada.

NOTA DE CONVENCIÓN (herramienta legado). Este validador y la evidencia JSON
del dossier de registro (2026-04-25) usan la convención de niveles
n = 1..N-1, bajo la cual la identidad global es

    mu_E + mu_A + mu_I = (N-1) * p^N.

La convención vigente del programa (Paper I y su bundle de reproducibilidad
en Zenodo, versión 1.3.1, concept DOI 10.5281/zenodo.20709159) es
n = 0..N-1: incluye el nivel raíz n = 0, cuya única bola se clasifica como
informativa, de modo que la identidad global pasa a ser

    mu_E + mu_A + mu_I = N * p^N,

con mu_E y mu_A sin cambio y mu_I aumentado exactamente en p^N. Ejemplo
(A. thaliana, N = 13, p = 2): mu_I pasa de 17648 a 25840 y la suma global
de 98304 a 106496 = 13 * 2^13.

Los valores esperados de este archivo se conservan a propósito en la
convención anterior porque validan la evidencia congelada del dossier. Para
los valores vigentes, la fuente de verdad es el bundle de Zenodo citado.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


@dataclass
class CaseSpec:
    name: str
    json_path: Path
    expected_scores: Dict[str, int]
    n: int
    p: int


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def check_case(case: CaseSpec) -> Tuple[bool, List[str]]:
    payload = load_json(case.json_path)
    scores = payload["scores"]
    issues: List[str] = []

    # Exact score checks
    for key, expected in case.expected_scores.items():
        got = scores.get(key)
        if got != expected:
            issues.append(f"{key}: esperado={expected}, obtenido={got}")

    # Identidad global en la convención legado (niveles n = 1..N-1):
    # mu_E + mu_A + mu_I = (N-1) * p^N. En la convención vigente (n = 0..N-1,
    # bundle Zenodo concept DOI 10.5281/zenodo.20709159) la identidad es
    # N * p^N. Ver la nota al inicio del archivo.
    left = scores["mu_E"] + scores["mu_A"] + scores["mu_I"]
    right = (case.n - 1) * (case.p**case.n)
    if left != right:
        issues.append(f"identidad_global: esperado={right}, obtenido={left}")

    # Identidad por nivel: E(n)+A(n)+I(n) = p^n, con n = 1..N-1 (convención
    # legado, sin el nivel raíz n = 0).
    E = scores.get("E", {})
    A = scores.get("A", {})
    I = scores.get("I", {})
    for level in range(1, case.n):
        k = str(level)
        got = int(E.get(k, 0)) + int(A.get(k, 0)) + int(I.get(k, 0))
        expected = case.p**level
        if got != expected:
            issues.append(f"identidad_nivel_{level}: esperado={expected}, obtenido={got}")

    return (len(issues) == 0), issues


def make_report(cases: List[CaseSpec], out_path: Path) -> None:
    lines = [
        "# Validación cruzada automática vs valores esperados",
        "",
        "Fecha: 2026-04-25",
        "",
        "Este reporte compara los JSON de evidencia con valores esperados de artículo y checks matemáticos.",
        "",
    ]

    global_ok = True
    for case in cases:
        ok, issues = check_case(case)
        status = "PASS" if ok else "FAIL"
        lines.append(f"## Caso: {case.name} — {status}")
        lines.append("")
        lines.append(f"- Archivo: `{case.json_path}`")
        lines.append(f"- Esperados: `{case.expected_scores}`")
        if ok:
            lines.append("- Resultado: todos los checks pasaron.")
        else:
            global_ok = False
            lines.append("- Issues:")
            for issue in issues:
                lines.append(f"  - {issue}")
        lines.append("")

    lines.append(f"## Resultado global: {'PASS' if global_ok else 'FAIL'}")
    lines.append("")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    root = Path(__file__).resolve()
    workspace = None
    for p in root.parents:
        if (p / "registro-dossier").exists():
            workspace = p
            break
    if workspace is None:
        raise RuntimeError("No se encontró raíz del workspace.")

    ev = workspace / "registro-dossier" / "06_evidencia_corridas"
    # Valores esperados en la convención legado n = 1..N-1, tal como quedaron
    # en la evidencia del dossier (2026-04-25). No actualizar a la convención
    # vigente: ver la nota al inicio del archivo.
    cases = [
        CaseSpec(
            name="Athaliana pi* (N=13, p=2)",
            json_path=ev / "core_scores_athaliana_pi_star.json",
            expected_scores={"mu_E": 26776, "mu_A": 53880, "mu_I": 17648},
            n=13,
            p=2,
        ),
        CaseSpec(
            name="C. elegans pi* (N=8, p=2)",
            json_path=ev / "core_scores_celegans_pi_star.json",
            expected_scores={"mu_E": 956, "mu_A": 392, "mu_I": 444},
            n=8,
            p=2,
        ),
        CaseSpec(
            name="Drosophila pi* (N=6, p=2)",
            json_path=ev / "core_scores_drosophila_pi_star.json",
            expected_scores={"mu_E": 64, "mu_A": 156, "mu_I": 100},
            n=6,
            p=2,
        ),
        CaseSpec(
            name="Toy identidad p=2, N=2",
            json_path=ev / "core_scores_p2.json",
            expected_scores={"mu_E": 0, "mu_A": 0, "mu_I": 4},
            n=2,
            p=2,
        ),
        CaseSpec(
            name="Toy identidad p=3, N=2",
            json_path=ev / "core_scores_p3.json",
            expected_scores={"mu_E": 0, "mu_A": 0, "mu_I": 9},
            n=2,
            p=3,
        ),
    ]
    out = workspace / "registro-dossier" / "06_evidencia_corridas" / "VALIDACION-CRUZADA-VALORES.md"
    make_report(cases, out)
    print(f"OK: reporte generado en {out}")


if __name__ == "__main__":
    main()
