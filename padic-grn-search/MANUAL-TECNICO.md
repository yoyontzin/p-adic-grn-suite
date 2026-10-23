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

### 3.1 Evaluación de la aptitud

Con `p = 2` y tabla completa, `MuCalculator` (núcleo) evalúa cada orden por
plegado de bits en `O(2^N)` operaciones de lista, con los mismos enteros que el
cálculo por bola de `compute_scores`. Tiempos medidos (Python 3.15, un núcleo):
0.7–0.9 ms por orden en *A. thaliana* (`N = 13`), unos 6 ms con `N = 16` y unos
33 ms con `N = 18`; el costo se duplica con cada gen adicional. Tras un
intercambio o una inserción solo se recalculan los niveles del bloque que cambia.
Con `p != 2` o tabla incompleta se usa el cálculo por bola.

### 3.2 Modo memético opcional

`pgsearch-ga --local-search` (o `random_ga_search(..., local_search=True)`)
aplica, a cada nuevo mejor orden de una generación, ascenso de colinas por
intercambios y movimientos de inserción, y lo devuelve a la población. No usa
números aleatorios, pero cambia la población; por eso no es el modo por
omisión. En *A. thaliana* (50 generaciones, población 30, semillas 1–20) el modo
por omisión llega a `mu_E = 26776` en 3 de 20 corridas y el memético en 20 de 20,
con el mismo tiempo (unos 0.65 s por corrida).

### 3.3 Redes grandes (`N >= 20`)

La memoria la domina la tabla de transiciones en Python: unos 100 MB medidos
con `N = 18`; por extrapolación (se cuadruplica cada dos genes), unos 400 MB con
`N = 20` y 1.6 GB con `N = 22`. Cada evaluación tarda, también por
extrapolación, unos 0.15 s con `N = 20` y 0.6 s con `N = 22`.
Estas corridas van en el servidor, nunca en la laptop:

```bash
rsync -a padic-grn-core padic-grn-search rtx1.cimat.mx:p-adic-grn-suite/
ssh rtx1.cimat.mx
tmux new -s ga
cd p-adic-grn-suite
python3 -m venv .venv && . .venv/bin/activate
pip install -e padic-grn-core -e padic-grn-search   # sin dependencias externas
pgsearch-ga red_N20.csv 20 2 --generations 200 --population 60 --seed 1 \
    --local-search --output-json ga_N20_seed1.json
# Ctrl+B D para dejarlo corriendo; tmux attach -t ga para revisarlo
```

Para `N >= 22` conviene el motor en JavaScript (`node`), o el programa dinámico
exacto mientras quepa en la memoria del servidor.

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
