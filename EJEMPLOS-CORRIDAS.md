# Ejemplos de corridas — p-adic GRN Suite

Este documento resume corridas representativas para la suite completa.

## Ejemplo 1: pipeline mínimo (toy p=2, N=2) con salida real

### Paso A. Scoring de ordering con core

```bash
cd software-suite/padic-grn-core
pgcore-score-ordering examples/transitions_p2_n2_identity.csv 2 2 "0 1" --output-json toy_core_scores.json
```

Archivo de salida:

- `toy_core_scores.json`

Salida real de referencia (`core_scores_p2.json`):

```json
{
  "mu_E": 0,
  "mu_A": 0,
  "mu_I": 4
}
```

### Paso B. Búsqueda con search

```bash
cd ../padic-grn-search
pgsearch-ga ../padic-grn-core/examples/transitions_p2_n2_identity.csv 2 2 --objective minimize_mu_E --generations 5 --population 10
```

Resultado real de referencia (`search_ga_p2.json`):

- `method`: `ga_pmx_tournament_elite`
- `best_ordering`: `[1, 0]`
- `best_fitness`: `0`
- `best_scores`: `mu_E=0`, `mu_A=0`, `mu_I=4`
- `evaluations`: `2`

### Paso C. Visualización con viz

```bash
cd ../padic-grn-viz
pgviz-dependency-graph
pgviz-hierarchical-graph
```

Resultado real actual:

- el módulo `viz` está funcional para generación de gráficas,
- en esta versión las capturas institucionales se reportan en:
  - `web-demo-institucional/screenshots/demo-ui-desktop.png`
  - `web-demo-institucional/screenshots/demo-ui-mobile.png`

## Ejemplo 2: validación completa de suite

```bash
cd software-suite
./RUN-ALL-VALIDATIONS.sh
```

Resultado real (última corrida reportada):

- ejecución de pruebas en los tres módulos,
- confirmación de estado global.

## Ejemplo 3: búsqueda exhaustiva (valores reales)

Resultados de referencia:

- `search_bbmin_p2.json`
  - `method`: `branch_and_bound_partial_muE`
  - `objective`: `minimize_mu_E`
  - `best_ordering`: `[0, 1]`
  - `best_fitness`: `0`
  - `best_scores`: `mu_E=0`, `mu_A=0`, `mu_I=4`
- `search_bbmax_p2.json`
  - `method`: `backtracking_exact`
  - `objective`: `maximize_mu_A`
  - `best_ordering`: `[0, 1]`
  - `best_fitness`: `0`
  - `best_scores`: `mu_E=0`, `mu_A=0`, `mu_I=4`

## Ejemplo 4: evidencias de referencia incluidas

Para revisión documental se incluyen salidas de referencia en:

- `registro-dossier/06_evidencia_corridas/core_scores_p2.json`
- `registro-dossier/06_evidencia_corridas/core_scores_p3.json`
- `registro-dossier/06_evidencia_corridas/search_ga_p2.json`
- `registro-dossier/06_evidencia_corridas/search_bbmin_p2.json`
- `registro-dossier/06_evidencia_corridas/search_bbmax_p2.json`

Estas salidas permiten verificar reproducibilidad contra corridas previas y se
usan como base documental en el PDF principal de la suite.
