#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
# Entorno unificado para imports entre paquetes de la suite.
export PYTHONPATH="$ROOT/padic-grn-core/src:$ROOT/padic-grn-search/src:$ROOT/padic-grn-viz/src:${PYTHONPATH:-}"

echo "[1/4] Tests core"
cd "$ROOT/padic-grn-core"
python3 -m pytest -q

echo "[2/4] Tests search"
cd "$ROOT/padic-grn-search"
python3 -m pytest -q

echo "[3/4] Tests viz"
cd "$ROOT/padic-grn-viz"
python3 -m pytest -q

echo "[4/4] Smoke web backend import"
cd "$ROOT/web-demo-institucional/backend"
python3 -c "import app; print('web_backend_import_ok')"

echo "Validaciones completadas."
