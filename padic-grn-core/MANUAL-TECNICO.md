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
3. Para cada nivel `n = 1, ..., N-1`, el espacio se particiona en cilindros/bolas `B_{n,m}`.
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

## 3. Arquitectura de implementación

- `src/padic_grn_core/cli.py`
  - `pgcore-score-ordering`
  - `pgcore-batch-scores`
- `src/padic_grn_core/adapters.py`
  - `from_grn_table(...)` (estable)
  - `from_discrete_system(...)` (interfaz extendida)

Motor de cálculo:

- La CLI se conecta al módulo validado de referencia:
  - `20250115-p-adic-GRN-rewrite/articulo-final/supplementary/code/core/compute_mu_generalized.py`

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

- Cálculo de una permutación: depende de `sum_{n=1}^{N-1} p^n` y del costo por bola.
- Para `N` grande, el cuello de botella está en evaluación masiva de bolas.
- El módulo soporta paralelismo en backend legado; la capa empaquetada lo expone indirectamente.

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
