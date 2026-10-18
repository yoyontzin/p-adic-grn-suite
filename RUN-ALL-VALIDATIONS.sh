#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
# Entorno unificado para imports entre paquetes de la suite.
export PYTHONPATH="$ROOT/padic-grn-core/src:$ROOT/padic-grn-search/src:$ROOT/padic-grn-viz/src:${PYTHONPATH:-}"

echo "[1/6] Tests core"
cd "$ROOT/padic-grn-core"
python3 -m pytest -q

echo "[2/6] Tests search"
cd "$ROOT/padic-grn-search"
python3 -m pytest -q

echo "[3/6] Tests viz"
cd "$ROOT/padic-grn-viz"
python3 -m pytest -q

echo "[4/6] Valores de referencia (niveles n = 0..N-1)"
cd "$ROOT/padic-grn-core"
python3 tools/validate_expected_values.py

echo "[5/6] Demos web contra fuerza bruta (requiere node)"
cd "$ROOT/web-demo-institucional"
python3 -m pytest -q tests

echo "[6/6] Smoke web backend import"
cd "$ROOT/web-demo-institucional/backend"
python3 -c "import app; print('web_backend_import_ok')"

echo "Validaciones completadas."
