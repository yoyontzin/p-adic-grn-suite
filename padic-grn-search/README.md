# padic-grn-search: Búsqueda y certificación de orderings p-ádicos

## Resumen

Paquete de optimización y certificación de ordenamientos (`GA` + `B&B`) sobre el modelo discreto definido en `padic-grn-core`.

## Introducción

`padic-grn-search` amplía el uso del núcleo matemático para explorar múltiples orderings de manera sistemática. Permite comparar candidatos y certificar soluciones en escenarios donde la búsqueda exhaustiva o heurística es necesaria.

## Comandos

- `pgsearch-ga`
- `pgsearch-bb-min-muE`
- `pgsearch-bb-max-muA`

## Implementación de optimización y certificación

La implementación actual incluye:

- `pgsearch-ga`: algoritmo genético con selección por torneo, cruce PMX,
  mutación por intercambio y élite.
- `pgsearch-bb-min-muE`: búsqueda por branch-and-bound/backtracking exacto para
  minimizar `mu_E`, con poda por cota parcial cuando está disponible.
- `pgsearch-bb-max-muA`: búsqueda exacta por backtracking para maximizar
  `mu_A`.

El cálculo de métricas se realiza con el núcleo matemático de referencia para
mantener continuidad y trazabilidad científica.

## Autoría oficial

1. Jesús Rogelio Pérez-Buendía (SECIHTI-CIMAT, Unidad Mérida) — ORCID `0000-0002-7739-4779`
2. Víctor Nopal-Coello (Universidad Autónoma del Estado de Hidalgo, UAEH) — ORCID `0000-0003-2608-3636`
3. Yuriria Cortés-Poza (IIMAS, UNAM — Unidad Académica Yucatán) — ORCID `0000-0003-1856-9206`

## DOI

- DOI del software: `DOI_SOFTWARE_SEARCH`.

## Financiamiento

Modelos matemáticos y computacionales no convencionales para el estudio y análisis de problemas relevantes en Biología, CF 2019/217367.
