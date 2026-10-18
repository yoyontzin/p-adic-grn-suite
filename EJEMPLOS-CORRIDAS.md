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

Salida real (`toy_core_scores.json`, campo `scores`):

```json
{
  "mu_E": 0,
  "mu_A": 0,
  "mu_I": 8
}
```

### Paso B. Búsqueda con search

```bash
cd ../padic-grn-search
pgsearch-ga ../padic-grn-core/examples/transitions_p2_n2_identity.csv 2 2 --objective minimize_mu_E --generations 5 --population 10
```

Resultado real:

- `method`: `ga_pmx_tournament_elite`
- `best_ordering`: `[1, 0]`
- `best_fitness`: `0`
- `best_scores`: `mu_E=0`, `mu_A=0`, `mu_I=8`
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
  - `best_scores`: `mu_E=0`, `mu_A=0`, `mu_I=8`
- `search_bbmax_p2.json`
  - `method`: `branch_and_bound_upper_muA`
  - `objective`: `maximize_mu_A`
  - `best_ordering`: `[0, 1]`
  - `best_fitness`: `0`
  - `best_scores`: `mu_E=0`, `mu_A=0`, `mu_I=8`

En el toy identidad cada nivel `n = 0, 1` es indiferente, así que
`mu_I = N p^N = 2 * 2^2 = 8`.

## Ejemplo 4: valores de referencia de las tres redes

```bash
cd padic-grn-core
python3 tools/validate_expected_values.py
```

Recalcula cada caso desde `examples/` y comprueba `mu_E + mu_A + mu_I = N p^N`
y `E(n) + A(n) + I(n) = p^n` para `n = 0..N-1`:

| Red | N | pi* | (mu_E, mu_A, mu_I) |
|---|---|---|---|
| *A. thaliana* | 13 | `[12,3,4,8,5,1,7,2,10,11,9,6,0]` | (26776, 53880, 25840) |
| *C. elegans* | 8 | `[1,0,2,3,6,5,4,7]` | (956, 392, 700) |
| *Drosophila* | 6 | `[0,5,3,4,2,1]` | (64, 156, 164) |

Las tablas de *A. thaliana* y *C. elegans* son las del bundle de
reproducibilidad de Paper II. La de *Drosophila* es la red de las demos web.
