# Manual técnico — padic-grn-core

## 1. Propósito científico

`padic-grn-core` implementa el núcleo matemático para cuantificar estabilidad jerárquica en sistemas dinámicos discretos finitos mediante métricas p-ádicas:

- `mu_E` (expansión),
- `mu_A` (contracción),
- `mu_I` (isometría/indiferencia),
- perfiles por nivel `E(n), A(n), I(n)`.

Este módulo es la base de `padic-grn-search` y `padic-grn-viz`.

## 2. Modelo matemático formal

Sea `f: C -> C` un sistema discreto determinista sobre `|C| = p^N`.

1. Se fija un ordenamiento `pi` de variables.
2. Se identifica cada estado con su expansión base `p`.
3. Para cada nivel `n = 0, ..., N-1`, el espacio se particiona en cilindros/bolas `B_{n,m}`.
   El nivel `n = 0` es una sola bola, todo el espacio; nunca es expansiva.
4. Para cada bola, se calcula `M_{n,m}` = primer nivel de discrepancia en imágenes bajo `f`.
5. Se clasifica:
   - `E` si `M_{n,m} < n`,
   - `I` si `M_{n,m} = n`,
   - `A` si `M_{n,m} > n`.
6. Con peso `w_n = p^(N-n)`:
   - `mu_E = sum_n E(n) * w_n`
   - `mu_I = sum_n I(n) * w_n`
   - `mu_A = sum_n A(n) * w_n`

Invariante de consistencia por nivel:

- `E(n) + I(n) + A(n) = p^n`.

Invariante global (cada nivel aporta `p^n * p^(N-n) = p^N`):

- `mu_E + mu_A + mu_I = N * p^N`.

## 3. Arquitectura de implementación

- `src/padic_grn_core/cli.py`
  - `pgcore-score-ordering`
  - `pgcore-batch-scores`
- `src/padic_grn_core/adapters.py`
  - `from_grn_table(...)` (estable)
  - `from_discrete_system(...)` (interfaz extendida)

Motor de cálculo:

- `src/padic_grn_core/compute_mu.py`, copia del motor del bundle de
  reproducibilidad de Paper II (`compute_mu_generalized.py`, commit `681aeee`).
  Lo usan la CLI y `padic-grn-search`.
- Cotas parciales para branch and bound, con prefijo `(a_0, ..., a_k)`:
  - `partial_mu_E(prefix, ...)`: suma de `E(n) p^(N-n)` para `n = 0..k`; cota inferior de `mu_E`.
  - `partial_mu_A(prefix, ...)`: suma de `A(n) p^(N-n)` para `n = 0..k`.
  - `mu_A_upper_bound(prefix, ...) = partial_mu_A + (N-1-k) p^N`; cota superior de `mu_A`.
- `tools/validate_expected_values.py` recalcula los valores de referencia desde
  `examples/` y comprueba ambas identidades.
- `tools/mu_dp_exact.js` (Node, p = 2) optimiza `mu_E` y `mu_A` sobre los `N!`
  ordenamientos por programación dinámica sobre subconjuntos.

## 4. Contrato de datos

Entrada mínima:

- `N` (entero positivo),
- `p` (entero >= 2),
- `transition_table.csv` con columnas:
  - `config_from`
  - `config_to`,
- `ordering` como permutación de `0..N-1`.

Salida estándar (`scores.json`):

- metadatos (`software`, `version`, `N`, `p`, `ordering`),
- `scores.mu_E`, `scores.mu_A`, `scores.mu_I`,
- perfiles por nivel.

## 5. Complejidad y límites prácticos

- Cálculo de una permutación: depende de `sum_{n=0}^{N-1} p^n` y del costo por bola.
- Para `N` grande, el cuello de botella está en evaluación masiva de bolas.
- `compute_scores(..., use_parallel=True)` reparte las bolas de cada nivel entre procesos.

## 6. Validación y trazabilidad

Validaciones obligatorias:

- ordering válido (permutación),
- CSV con cabeceras correctas,
- integridad de tabla de transición.

Trazabilidad:

- toda corrida debe guardar: parámetros, versión, ruta de datos y hash de entrada.

## 7. Extensibilidad técnica

1. Consolidar el motor en `src/padic_grn_core/engine.py`.
2. Añadir hash de entrada y checksum de salida.
3. Añadir pruebas deterministas para `p=2` y `p=3` en CI.
4. Añadir API python estable:
   - `compute_scores(...)`,
   - `validate_transition_table(...)`.
