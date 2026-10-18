# p-adic GRN Suite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19805643.svg)](https://doi.org/10.5281/zenodo.19805643)

Suite cientifica para analisis p-adico de redes de regulacion genetica en dinamica discreta finita.

## Componentes

- `padic-grn-core`: metricas p-adicas.
- `padic-grn-search`: optimizacion y certificacion de orderings.
- `padic-grn-viz`: visualizacion y reportes.
- `web-demo-institucional`: applets y evidencias visuales.

## Medidas y niveles

Para un ordenamiento de los `N` genes, cada nivel `n = 0, ..., N-1` parte el
espacio de `p^N` estados en `p^n` bolas, y cada bola se clasifica como
expansiva (E), atractiva (A) o indiferente (I). Con peso `p^(N-n)` por bola,

    mu_E + mu_A + mu_I = N * p^N,      E(n) + A(n) + I(n) = p^n.

El nivel `n = 0` es una sola bola (todo el espacio) y nunca es expansiva.
Valores de referencia: `padic-grn-core/tools/validate_expected_values.py`.
Historial de versiones: `CHANGELOG.md`.

## Cómo citar

- Versión actual: v1.2.1 (`CITATION.cff`).
- DOI conceptual (todas las versiones; resuelve a la más reciente):
  `10.5281/zenodo.19805643`, <https://doi.org/10.5281/zenodo.19805643>
- Releases: <https://github.com/yoyontzin/p-adic-grn-suite/releases>

La v1.0.1 (DOI de versión `10.5281/zenodo.19805644`, SWHID
`swh:1:dir:8441947af7175929272589c872a5dc63eb3c7fca`) es la versión depositada
del registro de autor y se conserva sin cambios.
