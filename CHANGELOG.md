# Changelog — p-adic GRN Suite

## v1.2.0 (2026-10-18)

### Medidas

- Todas las sumas por niveles recorren `n = 0, ..., N-1`, como en los artículos
  y en el bundle de reproducibilidad de Paper II: Python (`padic-grn-core`,
  `padic-grn-search`, `tools/`) y las dos demos web.
  - Identidades: `mu_E + mu_A + mu_I = N * p^N` y `E(n) + A(n) + I(n) = p^n`.
  - El nivel `n = 0` es una sola bola, todo el espacio, y nunca es expansiva,
    así que `mu_E` y los ordenamientos que lo minimizan no cambian respecto a
    v1.1.x. Esa bola suma `p^N` a `mu_A` si la raíz contrae y a `mu_I` si es
    isométrica.
  - *A. thaliana*, pi* = `[12,3,4,8,5,1,7,2,10,11,9,6,0]`: `(26776, 53880, 25840)`.
  - *C. elegans*, pi* = `[1,0,2,3,6,5,4,7]`: `(956, 392, 700)`.
  - *Drosophila*, pi* = `[0,5,3,4,2,1]`: `(64, 156, 164)`.

### padic-grn-core 0.2.0

- `src/padic_grn_core/compute_mu.py`: el motor queda dentro del paquete (copia
  de `compute_mu_generalized.py` del bundle de Paper II, commit `681aeee`). La
  CLI ya no busca el motor en carpetas externas al repositorio.
- Cotas parciales para branch and bound, sobre los niveles `n = 0..k` que fija
  un prefijo `(a_0, ..., a_k)`: `partial_mu_E`, `partial_mu_A` y
  `mu_A_upper_bound = partial_mu_A + (N-1-k) p^N`.
- `tools/validate_expected_values.py` recalcula los valores de referencia desde
  las tablas de `examples/`, en lugar de leer JSON externos.
- `examples/`: tablas de transición de *A. thaliana* y *C. elegans* (bundle de
  Paper II) y de *Drosophila* (la red de las demos web).
- `tools/mu_dp_exact.js`: programa dinámico exacto del bundle de Paper II.

### padic-grn-search 0.2.0

- `pgsearch-bb-max-muA` poda con la cota superior de `mu_A`; antes recorría el
  árbol completo.

### Demos web

- Perfiles, tablas, gráficas, slider y tiras de historia incluyen el nivel `n = 0`.
- `padic-grn-demo.html`: la tabla de *A. thaliana* tenía 6144 de los 8192 estados
  y no coincidía con la red publicada; ahora es la tabla completa del bundle
  (la misma de `padic-balls-animation.html`).
- Las bolas sin imágenes en la tabla cuentan como atractivas, igual que en el
  motor de Python, de modo que la identidad por nivel se mantiene.
