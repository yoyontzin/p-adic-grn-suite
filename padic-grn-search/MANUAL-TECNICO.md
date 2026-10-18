# Manual técnico — padic-grn-search

## 1. Propósito

`padic-grn-search` resuelve optimización y certificación de ordenamientos:

- optimización heurística (GA),
- certificación exacta/semiexacta (B&B, según configuración).

## 2. Arquitectura actual

- `src/padic_grn_search/engine.py`
  - `random_ga_search(...)` (GA con torneo + PMX + mutación + élite),
  - `branch_and_bound_search(...)` (certificación exacta por árbol de búsqueda),
  - `exhaustive_search(...)` (búsqueda completa, útil para verificación).
- `src/padic_grn_search/cli.py`
  - `pgsearch-ga`
  - `pgsearch-bb-min-muE`
  - `pgsearch-bb-max-muA`

La capa de búsqueda es nativa; el cálculo de métricas usa el núcleo matemático
de referencia para conservar continuidad científica.

## 3. Diseño algorítmico (resumen)

- GA nativo:
  - selección por torneo,
  - cruce PMX para permutaciones válidas,
  - mutación por intercambio,
  - élite por generación,
  - evaluación por objetivo (`mu_E` o `-mu_A`) con cache de evaluaciones.
- Branch-and-Bound / backtracking exacto:
  - para `minimize_mu_E` poda con la cota inferior `partial_mu_E` (niveles `n = 0..k` del prefijo),
  - para `maximize_mu_A` poda con la cota superior `partial_mu_A + (N-1-k) p^N`.
- Exhaustivo nativo:
  - recorre `N!` permutaciones,
  - opción de truncamiento por presupuesto de evaluaciones.

## 4. Entrada/Salida

Entrada:

- transiciones (`csv`),
- `N`, `p`,
- parámetros de corrida (`generations`, `population`, `seed`, `max-evaluations`).

Salida:

- mejor ordenamiento,
- métricas (`mu_E`, `mu_A`, `mu_I`),
- conteo de evaluaciones,
- conteo de ramas podadas (en B&B),
- indicador de truncamiento (si aplica).

## 5. Trazabilidad de ejecución

- fijar `--seed` en corridas de GA para reproducibilidad,
- registrar parámetros de entrada y objetivo,
- conservar JSON de salida de cada corrida como evidencia verificable.
