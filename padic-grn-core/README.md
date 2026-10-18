# padic-grn-core: Núcleo p-ádico para sistemas dinámicos discretos

## Resumen

Paquete base para calcular métricas p-ádicas (`mu_E`, `mu_A`, `mu_I`) en sistemas dinámicos discretos finitos.
Las medidas suman los niveles `n = 0, ..., N-1`, con `mu_E + mu_A + mu_I = N * p^N`.

## Introducción

`padic-grn-core` concentra el cálculo matemático central de la suite. Está orientado a redes de regulación genética multivaluadas y, mediante adaptadores, puede emplearse también en sistemas discretos generales.

## Comandos

- `pgcore-score-ordering`
- `pgcore-batch-scores`

## Uso rápido

```bash
pgcore-score-ordering transitions.csv 13 2 "12 4 3 8 1 5 7 2 10 11 9 6 0" --output-json scores.json
```

`pgcore-batch-scores` espera un CSV con columna `ordering`.

## Alcance objetivo

- Soporte oficial: GRN deterministas multivaluadas.
- Soporte experimental: sistemas discretos generales via adaptadores.

## Autoría oficial

1. Jesús Rogelio Pérez-Buendía (SECIHTI-CIMAT, Unidad Mérida) — ORCID `0000-0002-7739-4779`
2. Víctor Nopal-Coello (Universidad Autónoma del Estado de Hidalgo, UAEH) — ORCID `0000-0003-2608-3636`
3. Yuriria Cortés-Poza (IIMAS, UNAM — Unidad Académica Yucatán) — ORCID `0000-0003-1856-9206`

## DOI

- DOI del software: `DOI_SOFTWARE_CORE`.

## Financiamiento

Modelos matemáticos y computacionales no convencionales para el estudio y análisis de problemas relevantes en Biología, CF 2019/217367.
